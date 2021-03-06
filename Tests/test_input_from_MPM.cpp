#include <iostream>

#include "SandHypoplasticityByUmat.h"
#include "Tests_main.h"

void test_input_from_MPM()
{
	union
	{
		struct
		{
			size_t udstrain0[6];
			size_t udstrain1[6];
		};
		struct
		{
			double dstrain0[6];
			double dstrain1[6];
		};
	};

	const double *stress;
	MatModel::SandHypoplasticityByUmat shp;
	const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	const double R = 1.0e-4;
	const double ig_strain[6] = { -R / sqrt(3.0), -R / sqrt(3.0), -R / sqrt(3.0), 0.0, 0.0, 0.0 };
	shp.set_param(ini_stress, 0.817,
		33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		2.2, 1.1, R, 0.1, 5.5, ig_strain);
	shp.set_integration_step_ratio(0.05);

	udstrain0[0] = 13691036303143785760;
	udstrain0[1] = 13691036303143785760;
	udstrain0[2] = 13696328946189003372;
	udstrain0[3] = 0;
	udstrain0[4] = 0;
	udstrain0[5] = 0;
	
	//dstrain0[0] = 1.209773e-8;
	//dstrain0[1] = 1.209773e-8;
	//dstrain0[2] = -1.36460e-7;
	//dstrain0[3] = 4.345909e-27;
	//dstrain0[4] = -3.034195e-8;
	//dstrain0[5] = -3.034198e-8;

	udstrain1[0] = 13697637712029187180;
	udstrain1[1] = 13697637712029180616;
	udstrain1[2] = 13694206740447595366;
	udstrain1[3] = 13790794046982147559;
	udstrain1[4] = 13786290447354775423;
	udstrain1[5] = 4562918410500002898;

	const double *de0 = dstrain0;
	const double *de1 = dstrain1;

	//std::cout << dstrain0[0] << ", " << dstrain0[1] << ", "
	//		<< dstrain0[2] << ", " << dstrain0[3] << ", "
	//		<< dstrain0[4] << ", " << dstrain0[5] << "\n";
	shp.integrate(dstrain0);
	//std::cout << dstrain0[0] << ", " << dstrain0[1] << ", "
	//		<< dstrain0[2] << ", " << dstrain0[3] << ", "
	//		<< dstrain0[4] << ", " << dstrain0[5] << "\n";

	std::cout << "\n2nd integration\n";
	shp.integrate(dstrain1);
	stress = shp.get_stress();
	std::cout << "s: " << stress[0] << ", " << stress[1] << ", "
			  << stress[2] << ", " << stress[3] << ", "
			  << stress[4] << ", " << stress[5] << "\n";

	size_t end_tag = 0;
}
