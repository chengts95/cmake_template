#pragma once
#include "common.h"

class ICompManager
{

public:
	int indim = 0;
	virtual void init_dt(Real dt) = 0;
	virtual void update_Ihis(Real *current_res) = 0;
	virtual void update_Ieq(Real *prev_x, Real *current_res, Real t) = 0;
	virtual ICompManager *clone() const = 0;
	virtual void insert_MNA(Real *MNA, int ndim) = 0;
	virtual void insert_MNA_COO(COO& MNA) = 0;
	virtual void init_hisitems(int &total_his){};
	virtual void output_Ihis(Real *his){};
	virtual void scatter_Ihis(Real *his){};
	virtual ~ICompManager(){};
};

template <class _T>
class CompManager  : public ICompManager
{

	

public:
	std::vector<_T> m_Objects;
	CompManager (){};

	void add_comp(_T &&S)
	{

		m_Objects.emplace_back(S);
	}

	void add_comp(_T *S)
	{

		m_Objects.emplace_back(*S);
	}
	void insert_MNA(Real *MNA, int ndim)
	{
		indim = ndim;
		for (auto& i : m_Objects)
			i.insert_MNA(MNA, ndim);
	}
	void insert_MNA_COO(COO& MNA)
	{
		
		for (auto& i : m_Objects)
			i.insert_MNA(MNA, 1);
	}
	void init_dt(Real dt)
	{
		for (auto& i : m_Objects)
			i.init_dt(dt);
	}
	void update_Ihis(Real *current_res)
	{

		for (auto& i : m_Objects)
			i.update_Ihis(current_res);
	}
	void update_Ieq(Real *prev_x, Real *current_res, Real t)
	{

		for (int i = 0; i < m_Objects.size(); i++)
			m_Objects[i].update_Ieq(prev_x, current_res, t);
		// for (auto i : m_Objects)
		// 	i.update_Ieq(prev_x, current_res,t);
	}
};

