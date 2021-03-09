#include <fstream>

#include "SandHypoplasticityByUmat.h"
#include "Tests_main.h"

size_t integrate_time0 = 0;

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
	
	udstrain1[0] = 13697637712029187180;
	udstrain1[1] = 13697637712029180616;
	udstrain1[2] = 13694206740447595366;
	udstrain1[3] = 13790794046982147559;
	udstrain1[4] = 13786290447354775423;
	udstrain1[5] = 4562918410500002898;

	const double *de0 = dstrain0;
	const double *de1 = dstrain1;

	integrate_time0 = 1;
	shp.integrate(dstrain0);

	integrate_time0 = 2;
	shp.integrate(dstrain1);

	size_t end_tag = 0;
}
