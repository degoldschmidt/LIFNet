/****************************************************************************
 *	izhi.h																	*
 *																			*
 *  Created on:		Oct 04, 2014											*
 *  Author: 		Dennis Goldschmidt										*
 *  Email:			goldschmidtd@ini.phys.ethz.ch							*
 *  								 										*
 *  								 										*
 *  Copyright (C) 2014 by Dennis Goldschmidt								*
 *																			*
 *  This file is part of the program LIFNet									*
 *																			*
 *  LIFNet is free software: you can redistribute it and/or modify			*
 *  it under the terms of the GNU General Public License as published by	*
 *  the Free Software Foundation, either version 3 of the License, or		*
 *  (at your option) any later version.										*
 *																			*
 *  This program is distributed in the hope that it will be useful,			*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of			*
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the			*
 *  GNU General Public License for more details.							*
 *																			*
 *  You should have received a copy of the GNU General Public License		*
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.	*
 *																			*
 ***************************************************************************/

#ifndef IZHI_H_
#define IZHI_H_

#include <armadillo>
using namespace std;
using namespace arma;

class Izhi{
public:
	Izhi(int inNeurons){
		Ne = int(4*(inNeurons/5));
		Ni = int(inNeurons/5);
		N = Ne+Ni;

		input_train.zeros(N);
		external.zeros(N);
		noise.zeros(N);
		spikes.zeros(N);

		arma_rng::set_seed_random();
		re = randn(Ne);
		ri = randn(Ni);

		param_a = join_cols(0.02 * ones<vec>(Ne), 0.02*ones<vec>(Ni)+0.08*ri);
		param_b = join_cols(0.2*ones<vec>(Ne), 0.25*ones<vec>(Ni)-0.05*ri);
		param_c = join_cols(-65.*ones<vec>(Ne)+15.*(re%re), -65.*ones<vec>(Ni));
		param_d = join_cols(8.*ones<vec>(Ne)-6.*(re%re), 2.*ones<vec>(Ni));

		membr_potential = -65.*ones<vec>(N);
		membr_recover = param_b%membr_potential;

		W_rec = join_rows(0.5*(abs(randn(N,Ne))), -0.1*abs(randn(N,Ni)));
	};
	~Izhi(){
	};

	double I(int i){
		return input_train(i);
	}

	void set_input(double val, int index=0){
		external(index) = val;
	};


	void set_noiseE(double val){
		noise.subvec(0,Ne-1) = val*ones<vec>(Ne);
	};

	void set_noiseI(double val){
		noise.subvec(Ne,N-1) = val*ones<vec>(Ni);
	};

	vec Sp(){
		return spikes;
	};

	vec u(){
		return membr_recover;
	};

	double u(int i){
		return membr_recover(i);
	};

	vec du(){
		return param_a%(param_b%membr_potential-membr_recover);
	};

	double du(int i){
		return param_a(i)*(param_b(i)*membr_potential(i)-membr_recover(i));
	};

	void update(double dt){
		input_train = W_rec*spikes+external + noise;
		spikes.zeros(N);

		fired = find(membr_potential >= 30.);
		spikes(fired).ones();
		membr_potential(fired) = param_c(fired);
		membr_recover(fired)+= param_d(fired);

		membr_potential += dt*dV();
		membr_recover += dt*du();
	}

	vec V(){
		return membr_potential;
	};

	double V(int i){
		return membr_potential(i);
	};

	vec dV(){
		return 0.04*pow(membr_potential,2)+5.*membr_potential+140.-membr_recover+input_train;
	};

	double dV(int i){
		return 0.04*pow(membr_potential(i),2)+5.*membr_potential(i)+140.-membr_recover(i)+input_train(i);
	};

private:
	int N;
	int Ne;
	int Ni;
	mat W_rec;
	vec spikes;
	vec membr_potential;			// [mV]
	vec membr_recover;				// [mV]
	vec input_train;				// [A]
	vec external;
	vec noise;
	uvec fired;

	vec param_a;
	vec param_b;
	vec param_c;
	vec param_d;

	vec re;
	vec ri;
};



#endif /* IZHI_H_ */
