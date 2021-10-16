#include <iostream>
#include <vector>
#include <complex>
#include <cuComplex.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Функция освобождения памяти GPU (device_pointer_input — указатель на входное поле, device_pointer_output — на выходное поле (результат), device_pointer_x1 — на вектор x, device_pointer_x2 — на вектор y, device_pointer_x3 — на вектор u или z, device_pointer_x4 — на вектор v, device_pointer_parameters — на вектор параметров преобразования, device_pointer_dimension — на вектор размерностей, device_pointer_progress — на атомарную переменную прогресса процесса)
void freeGPUMemory(cuDoubleComplex* device_pointer_input, cuDoubleComplex* device_pointer_output, double* device_pointer_x1, double* device_pointer_x2, double* device_pointer_x3, double* device_pointer_x4, double* device_pointer_parameters, int* device_pointer_dimension, int* device_pointer_progress) {
    cudaFree(device_pointer_input);
    cudaFree(device_pointer_output);
    cudaFree(device_pointer_x1);
    cudaFree(device_pointer_x2);
    cudaFree(device_pointer_x3);
    cudaFree(device_pointer_x4);
    cudaFree(device_pointer_parameters);
    cudaFree(device_pointer_dimension);
    cudaFree(device_pointer_progress);
}

// Сложение двух CUDA double комплексных чисел (left — число слева, right — справа)
__device__ cuDoubleComplex operator+(const cuDoubleComplex& left, const cuDoubleComplex& right) {
    return cuCadd(left, right);
}

// Сложение двух CUDA double комплексных чисел с присвоением (left — переменная до знака "=", right — после)
__device__ cuDoubleComplex operator+=(cuDoubleComplex& left, const cuDoubleComplex& right) {
    left = left + right;
    return left;
}

// Перемножение двух CUDA double комплексных чисел (left — число слева, right — справа)
__device__ cuDoubleComplex operator*(const cuDoubleComplex& left, const cuDoubleComplex& right) {
    return cuCmul(left, right);
}

// Умножение CUDA double комплексного числа на вещественное double число (left — комплексное CUDA число, right — вещественное double число)
__device__ cuDoubleComplex operator*(const cuDoubleComplex& left, const double& right) {
    return cuCmul(left, make_cuDoubleComplex(right, 0));
}

// Умножение вещественного double числа на CUDA double комплексное число (left — вещественное double число, right — комплексное CUDA число)
__device__ cuDoubleComplex operator*(const double& left, const cuDoubleComplex& right) {
    return cuCmul(make_cuDoubleComplex(left, 0), right);
}

// Экспонента, принимающая в качестве аргумента CUDA double комплексное число
__device__ cuDoubleComplex exp(const cuDoubleComplex& value) {
    return exp(value.x) * make_cuDoubleComplex(cos(value.y), sin(value.y));
}

// Нахождение оптимального количества нитей в блоке (N — количество точек по одной оси)
int getNumberThreads(int N) {
    auto result = static_cast<int>(round(sqrt(N)));
    while (((N % result) != 0) || (result * result > 1024)) {
        result--;
    }
    return result;
}

// Функция отображения прогресса выполнения преобразования (now — сколько сейчас выполнено операций, max — общее кол-во операций)
__device__ void processing(int now, int max) {
    double percent;
    if (now == max) {
        percent = 100.;
    }
    else {
        percent = trunc(10000. * (static_cast<double>(now) / max)) / 100;
    }
    printf("\rВыполнено %2.2f%", percent);
}

// Обобщенное ядро вычисления интеграла Коллинза (output — массив результата (выходное поле), input — входное поле, x1 — вектор x, x2 — вектор y, x3 — вектор u или z, x4 — вектор v, parameters — массив параметров преобразования, dimension — массив размерностей, progress — атомарная переменная прогресса процесса вычисления интеграла Коллинза, transformType (вид преобразования: 0 — дробное преобразование Фурье, 1 — преобразование Френеля, 2 — другое (с заданной определенной ABCD матрицей)), OxyCrossSection (выбранное сечение: 1 — сечение в плоскости Ouv (поперечное), 0 — сечение в плоскости Ovz (продольное)))
__global__ void collinsKernel(cuDoubleComplex* output, const cuDoubleComplex* input, const double* x1, const double* x2, const double* x3, const double* x4, const double* parameters, const int* dimension, int* progress, int transformType, bool OxyCrossSection)
{
    auto pi = 3.14159265358979323846;
    auto q = blockIdx.x * blockDim.x + threadIdx.x;
    auto p = blockIdx.y * blockDim.y + threadIdx.y;
    auto hx = x1[1] - x1[0];
    auto hy = x2[0] - x2[1];
    auto wavelength = parameters[0];
    auto z = OxyCrossSection ? parameters[1] : x3[q];
    auto u = OxyCrossSection ? x3[q] : parameters[1];
    auto f = !transformType ? parameters[2] : 0;
    auto k = 2 * pi / wavelength;
    auto n1 = dimension[0];
    auto n2 = dimension[1];
    auto n3 = OxyCrossSection ? dimension[1] : dimension[2];
    auto A = 0.0;
    auto B = 0.0;
    auto D = 0.0;
    switch (transformType) {
    case 0:
        A = cos(pi * z / (2 * f));
        B = f * sin(pi * z / (2 * f));
        D = cos(pi * z / (2 * f));
        break;
    case 1:
        A = 1.;
        B = z;
        D = 1.;
        break;
    default:
        A = OxyCrossSection ? parameters[1] : parameters[2];
        B = OxyCrossSection ? parameters[2] : parameters[3];
        D = OxyCrossSection ? parameters[3] : parameters[4];
    }
    auto value = make_cuDoubleComplex(0, 0);
    for (auto i = 0; i < n1; i++) {
        for (auto j = 0; j < n1; j++) {
            auto arg = (k / (2 * B)) * (A * (x2[i] * x2[i] + x1[j] * x1[j]) - 2 * (x2[i] * x4[p] + x1[j] * u) + D * (x4[p] * x4[p] + u * u));
            value += input[i * n1 + j] * exp(make_cuDoubleComplex(0, arg));
        }
    }
    atomicAdd(progress, 1);
    processing(*progress, n2 * n3);
    output[p * n3 + q] = make_cuDoubleComplex(0, -(k / (2 * pi * B))) * value * hx * hy;
}

std::vector<std::vector<std::complex<double>>> calculateCollinsCUDA(const std::vector<std::vector<std::complex<double>>>& input, const std::vector<double>& x1, const std::vector<double>& x2, const std::vector<double>& x3, const std::vector<double>& x4, const std::vector<double>& parameters, const std::vector<int>& dimension, int transformType)
{
    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaError_t cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        throw std::runtime_error("cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
    }

    bool OxyCrossSection = dimension.size() == 3 ? false : true;

    // Allocate GPU buffers for vectors.

    auto n1 = dimension.at(0);
    auto n2 = dimension.at(1);
    auto n3 = OxyCrossSection ? dimension.at(1) : dimension.at(2);

    auto device_output = new cuDoubleComplex[n2 * n3];
    cuDoubleComplex* device_pointer_output = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_output, static_cast<unsigned long long>(n2) * n3 * sizeof(cuDoubleComplex));
    
    auto device_input = new cuDoubleComplex[input.size() * input.at(0).size()];
    for (auto i = 0; i < input.size(); i++) {
        for (auto j = 0; j < input.at(0).size(); j++) {
            device_input[i * n1 + j] = make_cuDoubleComplex(input.at(i).at(j).real(), input.at(i).at(j).imag());
        }
    }
    cuDoubleComplex* device_pointer_input = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_input, static_cast<unsigned long long>(n1) * n1 * sizeof(cuDoubleComplex));
    cudaStatus = cudaMemcpy(device_pointer_input, device_input, static_cast<unsigned long long>(n1) * n1 * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    
    std::vector<double> host_x1 = x1;
    double* pointer_x1 = host_x1.data();
    double* device_pointer_x1 = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_x1, n1 * sizeof(double));
    cudaStatus = cudaMemcpy(device_pointer_x1, pointer_x1, n1 * sizeof(double), cudaMemcpyHostToDevice);

    std::vector<double> host_x2 = x2;
    double* pointer_x2 = host_x2.data();
    double* device_pointer_x2 = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_x2, n1 * sizeof(double));
    cudaStatus = cudaMemcpy(device_pointer_x2, pointer_x2, n1 * sizeof(double), cudaMemcpyHostToDevice);

    std::vector<double> host_x3 = x3;
    double* pointer_x3 = host_x3.data();
    double* device_pointer_x3 = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_x3, n3 * sizeof(double));
    cudaStatus = cudaMemcpy(device_pointer_x3, pointer_x3, n3 * sizeof(double), cudaMemcpyHostToDevice);

    std::vector<double> host_x4 = x4;
    double* pointer_x4 = host_x4.data();
    double* device_pointer_x4 = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_x4, n2 * sizeof(double));
    cudaStatus = cudaMemcpy(device_pointer_x4, pointer_x4, n2 * sizeof(double), cudaMemcpyHostToDevice);

    std::vector<double> host_parameters = parameters;
    double* pointer_parameters = host_parameters.data();
    double* device_pointer_parameters = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_parameters, host_parameters.size() * sizeof(double));
    cudaStatus = cudaMemcpy(device_pointer_parameters, pointer_parameters, host_parameters.size() * sizeof(double), cudaMemcpyHostToDevice);

    std::vector<int> host_dimension = dimension;
    int* pointer_dimension = host_dimension.data();
    int* device_pointer_dimension = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_dimension, host_dimension.size() * sizeof(int));
    cudaStatus = cudaMemcpy(device_pointer_dimension, pointer_dimension, host_dimension.size() * sizeof(int), cudaMemcpyHostToDevice);
    
    int* device_pointer_progress = 0;
    cudaStatus = cudaMalloc((void**)&device_pointer_progress, sizeof(int));

    // Launch a kernel on the GPU with one thread for each element.
    dim3 threadsPerBlock(getNumberThreads(n3), getNumberThreads(n2));
    dim3 numBlocks(n3 / threadsPerBlock.x, n2 / threadsPerBlock.y);
    collinsKernel<<<numBlocks, threadsPerBlock>>>(device_pointer_output, device_pointer_input, device_pointer_x1, device_pointer_x2, device_pointer_x3, device_pointer_x4, device_pointer_parameters, device_pointer_dimension, device_pointer_progress, transformType, OxyCrossSection);
    
    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\ncollinsKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        freeGPUMemory(device_pointer_input, device_pointer_output, device_pointer_x1, device_pointer_x2, device_pointer_x3, device_pointer_x4, device_pointer_parameters, device_pointer_dimension, device_pointer_progress);
        throw std::runtime_error("Запуск ядра CUDA не удался!");
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "\ncudaDeviceSynchronize returned error code %d after launching collinsKernel!\n", cudaStatus);
        freeGPUMemory(device_pointer_input, device_pointer_output, device_pointer_x1, device_pointer_x2, device_pointer_x3, device_pointer_x4, device_pointer_parameters, device_pointer_dimension, device_pointer_progress);
        throw std::runtime_error("Синхронизация данных между хостом и устройством завершилась неудачей!");
    }

    std::cout << "\rВыполнено 100.00%" << std::endl;

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(device_output, device_pointer_output, static_cast<unsigned long long>(n2) * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        freeGPUMemory(device_pointer_input, device_pointer_output, device_pointer_x1, device_pointer_x2, device_pointer_x3, device_pointer_x4, device_pointer_parameters, device_pointer_dimension, device_pointer_progress);
        throw std::runtime_error("Копирование результата в ОЗУ завершилось неудачей!");
    }

    freeGPUMemory(device_pointer_input, device_pointer_output, device_pointer_x1, device_pointer_x2, device_pointer_x3, device_pointer_x4, device_pointer_parameters, device_pointer_dimension, device_pointer_progress);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        throw std::runtime_error("cudaDeviceReset failed!");
    }

    std::vector<std::vector<std::complex<double>>> result;
    result.reserve(x4.size());
    for (auto i = 0; i < x4.size(); i++) {
        auto row = std::vector<std::complex<double>>();
        row.reserve(x4.size());
        for (auto j = 0; j < x3.size(); j++) {
            row.emplace_back(std::complex<double>(device_output[i * n3 + j].x, device_output[i * n3 + j].y));
        }
        result.emplace_back(row);
    }

    return result;
}