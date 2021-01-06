#include <cstdlib>
#include "MaterialModelInitializer.h"
#include "Tests_main.h"

int main(int argc, char **argv)
{
	MatModel::MaterialModelInitializer::init();

	test_SandHypoplasticityByUmat();

	system("pause");
	return 0;
}
