#include <fstream>

#include "SandHypoplasticityByUmat.h"
#include "Tests_main.h"

void test_SandHypoplasticityByUmat()
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		Consolidation = 2
	};

	AnalysisType tp = AnalysisType::TriaxialUndrained;
	double de = -0.1;
	size_t inc_num = 1000;
	double ini_stress[6] = { -1.2e6, -1.2e6, -1.2e6, 0.0, 0.0, 0.0 };
	//double ini_stress[6] = { -0.9e6, -0.9e6, -0.9e6, 0.0, 0.0, 0.0 };
	double e0 = 0.67; // loose
	//double e0 = 0.53; // dense

	std::fstream res_file;
	const double(*Dep_mat)[6];
	double dstrain[6];
	const double *stress;

	de = de / 1000;
	
	MatModel::SandHypoplasticityByUmat shp;
	shp.set_param(ini_stress, e0,
		33.0, 1.6e9, 0.19, 0.44, 0.85, 1.00, 0.25, 1.0,
		5.0, 2.0, 1.0e-4, 0.5, 6.0);
	res_file.open("shp_res.csv", std::ios::out | std::ios::binary);
	stress = shp.get_stress();
	res_file << "strain, s11, s22, s33, q, p, e\n"
		<< 0.0 << ", "
		<< stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
		<< stress[0] - stress[1] << ", "
		<< (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
		<< shp.get_e() << "\n";

	for (size_t i = 0; i < inc_num; ++i)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			Dep_mat = reinterpret_cast<const double(*)[6]>(shp.get_Dep_mat());
			dstrain[0] = de;
			dstrain[1] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de;
			dstrain[2] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = de;
			dstrain[1] = -0.5 * de;
			dstrain[2] = -0.5 * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::Consolidation:
			dstrain[0] = de;
			dstrain[1] = 0.0;
			dstrain[2] = 0.0;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		default:
			res_file.close();
			return;
		}

		int res = shp.integrate(dstrain);

		// output result
		stress = shp.get_stress();
		res_file << de * (i + 1) << ", "
				 << stress[0] << ", "
				 << stress[1] << ", "
				 << stress[2] << ", "
				 <<  stress[0] - stress[1] << ", "
				 << (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
				 << shp.get_e() << "\n";
	}

	res_file.close();
}
