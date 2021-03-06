#ifndef __Material_Model_h__
#define __Material_Model_h__

#include "LinkList.hpp"

namespace MatModel
{
	class MaterialModel;
	template <typename MModel> class __Mat_Model_Container__;

	typedef int(*CMIntFunc)(MaterialModel *_self, double dstrain[6]);
	
	class MaterialModel
	{
	protected:
		const char *type;
		size_t id;

		union // stress
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		union // stress increment
		{
			double dstress[6];
			struct { double ds11, ds22, ds33, ds12, ds23, ds31; };
		};
		union // elastic strain increment
		{
			double dstrain_e[6];
			struct { double dee11, dee22, dee33, dee12, dee23, dee31; };
		};
		union // plastic strain increment
		{
			double dstrain_p[6];
			struct { double dep11, dep22, dep33, dep12, dep23, dep31; };
		};
		union // Elastic stiffness
		{
			double De_mat[6][6];
			double De_mat_array[36];
		};
		union // Elasto-plastic stiffness
		{
			double Dep_mat[6][6];
			double Dep_mat_array[36];
		};

		CMIntFunc integration_func;

	public:
		MaterialModel(
			CMIntFunc _inte_func = nullptr,
			const char *_type = "MaterialModel"
			) :
			integration_func(_inte_func), type(_type) {}
		~MaterialModel() {}

		inline int integrate(double dstrain[6]) { return (*integration_func)(this, dstrain); }

		inline size_t get_id() const { return id; }
		inline void set_id(size_t _id) { id = _id; }

		inline const char *get_type() { return type; }
		inline const double *get_stress()    { return stress; }
		inline const double *get_dstress()   { return dstress; }
		inline const double *get_dstrain_e() { return dstrain_e; }
		inline const double *get_dstrain_p() { return dstrain_p; }
		inline const double *get_De_mat()  { return De_mat_array; }
		inline const double *get_Dep_mat() { return Dep_mat_array; }

		union
		{
			void *ext_data_pt;
			unsigned long long ext_data_ull;
			long long ext_data_ll;
			double ext_data_d;
		};

	protected:
		template <typename MModel> friend class __Mat_Model_Container__;
		LinkListPointer<MaterialModel> pointer_by_container;
	};
};

#endif