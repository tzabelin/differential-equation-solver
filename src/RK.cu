#include "cuda_runtime.h"
#include "device_launch_parameters.h"
namespace Runge_Kutta
{
template<typename System, typename precision>
__device__ void Order_2(System system, precision *args, const precision dx, const int n, const int args_number)
{
	precision *args_tmp = new precision[args_number];

	// k is different derevatives in different points for calculation 
	precision *k1 = new precision[args_number - 1];
	precision *k2 = new precision[args_number - 1];

	system(args, k1);

	args_tmp[0] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx * k1[i - 1];
	}
	system(args_tmp, k2);

	args[n*args_number] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args[n*args_number + i] = args[i] + dx / 2 * k1[i - 1] + dx / 2 * k2[i - 1];
	}

	delete[] args_tmp;
	delete[] k1;
	delete[] k2;
}

template<typename System, typename precision>
__device__ void Order_3(System system, precision *args, const precision dx, const int n, const int args_number)
{
	precision *args_tmp = new precision[args_number];

	// k is different derevatives in different points for calculation
	precision *k1 = new precision[args_number - 1];
	precision *k2 = new precision[args_number - 1];
	precision *k3 = new precision[args_number - 1];

	system(args, k1);

	args_tmp[0] = args[0] + dx / 2;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx / 2 * k1[i - 1];
	}
	system(args_tmp, k2);

	args_tmp[0] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] - dx * k1[i - 1] + 2 * dx * (k2[i - 1]);
	}
	system(args_tmp, k3);

	args[n*args_number] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args[n*args_number + i] = args[i] + dx / 6 * k1[i - 1] + dx * 2 / 3 * k2[i - 1] + dx / 6 * k3[i - 1];
	}

	delete[] args_tmp;
	delete[] k1;
	delete[] k2;
	delete[] k3;
}

template<typename System, typename precision>
__device__ void Order_4(System system, precision *args, const precision dx, const int n, const int args_number)
{
	precision *args_tmp = new precision[args_number];

	// k is different derevatives in different points for calculation
	precision *k1 = new precision[args_number - 1];
	precision *k2 = new precision[args_number - 1];
	precision *k3 = new precision[args_number - 1];
	precision *k4 = new precision[args_number - 1];

	system(args, k1);

	args_tmp[0] = args[0] + dx / 2;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx / 2 * k1[i - 1];
	}
	system(args_tmp, k2);

	args_tmp[0] = args[0] + dx / 2;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx / 2 * (k2[i - 1]);
	}
	system(args_tmp, k3);

	args_tmp[0] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx * (k3[i - 1]);
	}
	system(args_tmp, k4);

	args[n*args_number] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args[n*args_number + i] = args[i] + dx / 6 * k1[i - 1] + dx / 3 * k2[i - 1] + dx / 3 * k3[i - 1] + dx / 6 * k4[i - 1];
	}

	delete[] args_tmp;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
}

template<typename System, typename precision>
__device__ void Order_8(System system, precision* args, const precision dx, const int n, const int args_number)
{
	precision* args_tmp = new precision[args_number];
	precision* k = new precision[args_number * 13]; 
	const precision a[13][12] = {
	{0},
	{1.0 / 18.0},
	{1.0 / 48.0, 1.0 / 16.0},
	{1.0 / 32.0, 0, 3.0 / 32.0},
	{5.0 / 16.0, 0, -75.0 / 64.0, 75.0 / 64.0},
	{3.0 / 80.0, 0, 0, 3.0 / 16.0, 3.0 / 20.0},
	{29443841.0 / 614563906.0, 0, 0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0},
	{16016141.0 / 946692911.0, 0, 0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0},
	{39632708.0 / 573591083.0, 0, 0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0},
	{246121993.0 / 1340847787.0, 0, 0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0},
	{-1028468189.0 / 846180014.0, 0, 0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0},
	{185892177.0 / 718116043.0, 0, 0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0},
	{403863854.0 / 491063109.0, 0, 0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0, 0}
	};
	const precision b[13] = { 14005451.0 / 335480064.0, 0, 0, 0, 0, -59238493.0 / 1068277825.0, 181606767.0 / 758867731.0, 561292985.0 / 797845732.0, -1041891430.0 / 1371343529.0, 760417239.0 / 1151165299.0, 118820643.0 / 751138087.0, -528747749.0 / 2220607170.0, 1.0 / 4.0 };
	const precision c[13] = { 0.0, 1.0 / 18.0, 1.0 / 12.0, 1.0 / 8.0, 5.0 / 16.0, 3.0 / 8.0, 59.0 / 400.0, 93.0 / 200.0, 5490023248.0 / 9719169821.0, 13.0 / 20.0, 1201146811.0 / 1299019798.0, 1.0, 1.0 };

	system(args, k);

	for (int j = 0; j < 13; j++) {
		args_tmp[0] = args[0] + c[j] * dx; 
		for (int i = 1; i < args_number; i++) {
			precision sum = 0;
			for (int l = 0; l < j; l++) {
				sum += a[j][l] * k[l * args_number + i - 1];
			}
			args_tmp[i] = args[i] + dx * sum;
		}
		system(args_tmp, &k[j * args_number]);
	}

	args[n * args_number] = args[0] + dx;
	for (int i = 1; i < args_number; i++) {
		precision sum = 0;
		for (int j = 0; j < 13; j++) {
			sum += b[j] * k[j * args_number + i - 1];
		}
		args[n * args_number + i] = args[i] + dx * sum;
	}

	delete[] args_tmp;
	delete[] k;
}
}
namespace Adams_Bashforth {
	template<typename System, typename precision>
	__device__ void Adams_Bashforth_2(System system, precision* args, precision* previous_k, const precision dx, const int n, const int args_number) {
		precision* k = new precision[args_number - 1];
		precision* args_tmp = new precision[args_number];

		// Compute current derivatives
		system(args, k);

		// Update the arguments based on the Adams-Bashforth formula
		args[n * args_number] = args[0] + dx;
		for (int i = 1; i < args_number; i++) {
			args[n * args_number + i] = args[i] + dx * ((3.0 / 2.0) * k[i - 1] - (1.0 / 2.0) * previous_k[i - 1]);
		}

		// Update previous derivatives for next step
		for (int i = 0; i < args_number - 1; i++) {
			previous_k[n * (args_number-1) + i] = k[i];
		}


		delete[] k;
	}

	//To make the first step before starting Adams-Bashforth
	template<typename System, typename precision>
	__device__ void RK_8(System system, precision* args, precision* previous_k, const precision dx, const int n, const int args_number)
	{
		precision* args_tmp = new precision[args_number];
		precision* k = new precision[args_number * 13];
		const precision a[13][12] = {
		{0},
		{1.0 / 18.0},
		{1.0 / 48.0, 1.0 / 16.0},
		{1.0 / 32.0, 0, 3.0 / 32.0},
		{5.0 / 16.0, 0, -75.0 / 64.0, 75.0 / 64.0},
		{3.0 / 80.0, 0, 0, 3.0 / 16.0, 3.0 / 20.0},
		{29443841.0 / 614563906.0, 0, 0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0},
		{16016141.0 / 946692911.0, 0, 0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0},
		{39632708.0 / 573591083.0, 0, 0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0},
		{246121993.0 / 1340847787.0, 0, 0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0},
		{-1028468189.0 / 846180014.0, 0, 0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0},
		{185892177.0 / 718116043.0, 0, 0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0},
		{403863854.0 / 491063109.0, 0, 0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0, 0}
		};
		const precision b[13] = { 14005451.0 / 335480064.0, 0, 0, 0, 0, -59238493.0 / 1068277825.0, 181606767.0 / 758867731.0, 561292985.0 / 797845732.0, -1041891430.0 / 1371343529.0, 760417239.0 / 1151165299.0, 118820643.0 / 751138087.0, -528747749.0 / 2220607170.0, 1.0 / 4.0 };
		const precision c[13] = { 0.0, 1.0 / 18.0, 1.0 / 12.0, 1.0 / 8.0, 5.0 / 16.0, 3.0 / 8.0, 59.0 / 400.0, 93.0 / 200.0, 5490023248.0 / 9719169821.0, 13.0 / 20.0, 1201146811.0 / 1299019798.0, 1.0, 1.0 };

		system(args, k);

		for (int j = 0; j < 13; j++) {
			args_tmp[0] = args[0] + c[j] * dx;
			for (int i = 1; i < args_number; i++) {
				precision sum = 0;
				for (int l = 0; l < j; l++) {
					sum += a[j][l] * k[l * args_number + i - 1];
				}
				args_tmp[i] = args[i] + dx * sum;
			}
			system(args_tmp, &k[j * args_number]);
		}

		args[n * args_number] = args[0] + dx;
		for (int i = 1; i < args_number; i++) {
			precision sum = 0;
			for (int j = 0; j < 13; j++) {
				sum += b[j] * k[j * args_number + i - 1];
			}
			args[n * args_number + i] = args[i] + dx * sum;
		}
		// Update previous derivatives for next step
		for (int i = 0; i < args_number - 1; i++) {
			previous_k[n * (args_number - 1)] = k[i];
		}
		delete[] args_tmp;
		delete[] k;
	}
}