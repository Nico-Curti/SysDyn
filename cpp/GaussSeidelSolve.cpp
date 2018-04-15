#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <cstring>

float* GaussSeidel(float **A, float *b, const int &N, const int &M)
{
	const int ITERATION_LIMIT = (int)1e20;
	int cnt;
	float *x_new = new float[M],
		  *res = new float[M],
			diff;
	for(int it = 0; it < ITERATION_LIMIT; ++it)
	{
#ifdef DEBUG
		std::cout << "Current solution: ";
		std::copy(x_new, x_new + M, std::ostream_iterator<float>(std::cout, " "));
		std::cout << std::endl;
#endif
		cnt = 0;
		for(int i = 0; i < N; ++i)
		{
			diff = std::inner_product(A[i], A[i] + i, x_new, 0.f) + 
				   std::inner_product(A[i] + i + 1, A[i] + M, x_new + i + 1, 0.f);
			x_new[i] = (b[i] - diff) / A[i][i];
			cnt += (std::fabs(res[i] - x_new[i]) > 1e-8) ? 1 : 0;
		}
		if(!cnt) break;
		std::memcpy(x_new, res, sizeof(float)*M);
	}
	delete[] x_new;
	return res;
}


int main(int argc, char **argv)
{
	float **A = new float*[4],
		   *b = new float[4],
		   *res = new float[4],
		   *errors = new float[4];
	std::generate(A, A + 4, [](){return new float[4];});
	A[0][0] = 10.f; A[0][1] = -1.f; A[0][2] = 2.f;  A[0][3] = 0.f;
    A[1][0] = -1.f; A[1][1] = 11.f; A[1][2] = -1.f; A[1][3] = 3.f;
    A[2][0] = 2.f;  A[2][1] = -1.f; A[2][2] = 10.f; A[2][3] = -1.f;
    A[3][0] = 0.f;  A[3][1] = 3.f;  A[3][2] = -1.f; A[3][3] = 8.f;

    b[0] = 6.f; b[1] = 25.f; b[2] = -11.f; b[3] = 15.f;

    std::cout << "System:" << std::endl;
    for(int i = 0; i < 4; ++i)
    {
    	for(int j = 0; j < 3; ++j) std::cout << A[i][j] << "*x" << j << " + ";
    	std::cout << A[i][3] << "*x" << 3 << " = " << b[i] << std::endl;
    }

    res = GaussSeidel(A, b, 4, 4);
    std::cout << "Solution:" << std::endl;
    std::copy(res, res + 4, std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;
    std::transform(A, A + 4, b, errors, [&res](const float *Ai, const float &bi)
    								 {
    								 	return std::inner_product(Ai, Ai + 4, res, 0.f) - bi;
    								 });
    std::cout << "Errors:" << std::endl;
    std::copy(errors, errors + 4, std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;

    for(int i = 0; i < 4; ++i) delete[] A[i];
    delete[] A;
	delete[] b;
	delete[] res;
	delete[] errors;

	return 0;
}
