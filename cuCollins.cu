#include <iostream>
#include <vector>
#include <complex>
#include "CollinsIntegralGPU.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuComplex.h>

using namespace std;

cudaError_t collinsWithCuda(cuDoubleComplex* input, cuDoubleComplex* output, double* x, double* y, double* u, double* v, int n1, int n2, double* fieldParameters);

__global__ void collinsKernel(cuDoubleComplex* output, cuDoubleComplex* input, double* x, double* y, double* u, double* v, int n1, int n2, double hx, double hy, double A, double B, double D, double k)
{
    int q = blockIdx.x * blockDim.x + threadIdx.x;
    int p = blockIdx.y * blockDim.y + threadIdx.y;
    cuDoubleComplex value = make_cuDoubleComplex(0, 0);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            double arg = (k / (2 * B)) * (A * (y[i] * y[i] + x[j] * x[j]) - 2 * (y[i] * v[p] + x[j] * u[q]) + D * (v[p] * v[p] + u[q] * u[q]));
            value = cuCadd(value, cuCmul(input[i * n1 + j], make_cuDoubleComplex(cos(arg), sin(arg))));
        }
    }
    output[p * n2 + q] = cuCmul(make_cuDoubleComplex(0, -(k / (2 * 3.14159265358979323846 * B))), cuCmul(value, make_cuDoubleComplex(hx * hy, 0)));
}

vector<vector<complex<double>>> calculateCollinsCUDA(vector<vector<complex<double>>>& inputFunction, vector<double>& x1, vector<double>& x2, vector<double>& x3, vector<double>& x4, int n1, int n2, double waveNumber, vector<double> limits, vector<vector<double>> matrixABCD)
{
    auto input = new cuDoubleComplex[inputFunction.size() * inputFunction.at(0).size()];
    for (auto i = 0; i < inputFunction.size(); i++) {
    	for (auto j = 0; j < inputFunction.at(0).size(); j++) {
    		input[i * n1 + j] = make_cuDoubleComplex(inputFunction.at(i).at(j).real(), inputFunction.at(i).at(j).imag());
    	}
    }

    auto output = new cuDoubleComplex[inputFunction.size() * inputFunction.at(0).size()];

    auto x = new double[x1.size()];
	for (auto i = 0; i < x1.size(); i++) {
		x[i] = x1.at(i);
	}

	auto y = new double[x2.size()];
	for (auto i = 0; i < x2.size(); i++) {
		y[i] = x2.at(i);
	}

	auto u = new double[x3.size()];
	for (auto i = 0; i < x3.size(); i++) {
		u[i] = x3.at(i);
	}

	auto v = new double[x4.size()];
	for (auto i = 0; i < x4.size(); i++) {
		v[i] = x4.at(i);
	}

    auto fieldParameters = new double[6];
    fieldParameters[0] = 2 * limits.at(0) / n1;
    fieldParameters[1] = 2 * limits.at(1) / n1;
    fieldParameters[2] = matrixABCD.at(0).at(0);
    fieldParameters[3] = matrixABCD.at(0).at(1);
    fieldParameters[4] = matrixABCD.at(1).at(1);
    fieldParameters[5] = waveNumber;

    // Add vectors in parallel.
    cudaError_t cudaStatus = collinsWithCuda(input, output, x, y, u, v, n1, n2, fieldParameters);
    if (cudaStatus != cudaSuccess) {
        error("collinsWithCuda failed!");
    }

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        error("cudaDeviceReset failed!");
    }

    vector<vector<complex<double>>> result;
    for (auto i = 0; i < inputFunction.size(); i++) {
        result.push_back(vector<complex<double>>());
        for (auto j = 0; j < inputFunction.at(0).size(); j++) {
            result.back().push_back(complex<double>(output[i * n2 + j].x, output[i * n2 + j].y));
        }
    }

    return result;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t collinsWithCuda(cuDoubleComplex* input, cuDoubleComplex* output, double* x, double* y, double* u, double* v, int n1, int n2, double* fieldParameters)
{
    cuDoubleComplex* dev_in = 0;
    cuDoubleComplex* dev_out = 0;
    size_t pitch_in;
    size_t pitch_out;
    double* dev_x = 0;
    double* dev_y = 0;
    double* dev_u = 0;
    double* dev_v = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output).
    cudaStatus = cudaMallocPitch((void**)&dev_in, &pitch_in, n1 * sizeof(cuDoubleComplex), n1);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMallocPitch((void**)&dev_out, &pitch_out, n2 * sizeof(cuDoubleComplex), n2);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_x, n1 * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_y, n1 * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_u, n2 * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_v, n2 * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy2D(dev_in, pitch_in, input, n1 * sizeof(cuDoubleComplex), n1 * sizeof(cuDoubleComplex), n1, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_x, x, n1 * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_y, y, n1 * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_u, u, n2 * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_v, v, n2 * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks(n2 / threadsPerBlock.x, n2 / threadsPerBlock.y);
    collinsKernel << <numBlocks, threadsPerBlock >> > (dev_out, dev_in, dev_x, dev_y, dev_u, dev_v, n1, n2, fieldParameters[0], fieldParameters[1], fieldParameters[2], fieldParameters[3], fieldParameters[4], fieldParameters[5]);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(output, dev_out, n2 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_in);
    cudaFree(dev_out);
    cudaFree(dev_x);
    cudaFree(dev_y);
    cudaFree(dev_u);
    cudaFree(dev_v);

    return cudaStatus;
}