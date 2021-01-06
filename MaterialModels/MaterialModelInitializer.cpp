#include "MaterialModels_pcp.h"

#include "SandHypoplasticityByUmat.h"
#include "SandHypoplasticityCavitationByUmat.h"
#include "MaterialModelInitializer.h"

namespace MatModel
{

MaterialModelInitializer::MaterialModelInitializer() {}

MaterialModelInitializer::~MaterialModelInitializer()
{
	if (is_init)
	{
		SandHypoplasticityByUmat::free();
		SandHypoplasticityCavitationByUmat::free();
		is_init = false;
	}
}

bool MaterialModelInitializer::is_init = false;

MaterialModelInitializer MaterialModelInitializer::instance;

int MaterialModelInitializer::init()
{
	if (is_init) return 0;

	is_init = true;
	SandHypoplasticityByUmat::init(
		"../SandHypoplasticity/Debug/SandHypoplasticity.dll");
	SandHypoplasticityCavitationByUmat::init(
		"../SandHypoplasticityCavitation/Debug/SandHypoplasticityCavitation.dll");
	return 0;
}

}