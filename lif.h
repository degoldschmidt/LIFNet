/*****************************************************************************
 *  lif.h                                                                    *
 *                                                                           *
 *  Created on:   Oct 04, 2014                                               *
 *  Author:       Dennis Goldschmidt                                         *
 *  Email:        goldschmidtd@ini.phys.ethz.ch                              *
 *                                                                           *
 *                                                                           *
 *  Copyright (C) 2014 by Dennis Goldschmidt								 *
 *  																		 *
 *  This file is part of the program NaviSim                                 *
 *                                                                           *
 *  NaviSim is free software: you can redistribute it and/or modify     	 *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.    *
 *                                                                           *
 ****************************************************************************/

#ifndef LIF_H_
#define LIF_H_

#include <armadillo>
using namespace std;
using namespace arma;

class LIF {
public:
	LIF(int inNeurons, double inResis = 10., double inCapac=1., double inRefra=4.){
		N = inNeurons;
		membr_potential = V_rest * ones<vec>(N);
		input_train.zeros(N);
		resistance = inResis * ones<vec>(N);
		capacitance = inCapac * ones<vec>(N);
		connectivity_mat.zeros(N, N);
		t_const = inResis * inCapac;
		t_reset.zeros(N);
		external.zeros(N);
		noise.zeros(N);
		refract_period = inRefra;
		fired.zeros(N);
	};
	~LIF(){

	};

	double Cm(int i=0){
		double out=as_scalar(capacitance(i));
		return out;};
	double I(int i=0){
		double out=as_scalar(input_train(i));
		return out;};
	double Rm(int i=0){
		double out=as_scalar(resistance(i));
		return out;};
	void set_input(double val, int index=0){
		external(index) = val;
	};
	void set_noise(double val){
		noise = val*ones<vec>(N);
	};

	void update(double t, double dt){
		input_train = /*connectivity_mat * membr_potential +*/ external + noise;
		for(int i=0; i < fired.n_elem; i++){
			if(fired(i) == 1){
				membr_potential(ind(i)) = V_rest;
				fired(i) = 0;
			}
			else
				membr_potential =  membr_potential + (dt/t_const)*(V_rest*ones<vec>(N) - membr_potential + input_train%resistance);
		}

		vec thr_vec = spike_thres*ones<vec>(N);
		ind = find(membr_potential >= thr_vec);
		for(int i=0; i < ind.n_elem; i++){
			if(t > t_reset(ind(i))){
				membr_potential(ind(i)) = spike;
				t_reset(ind(i)) = t + refract_period;
				fired(ind(i)) = 1;
			}
		}
	};
	double Vm(int i=0){
		double out=as_scalar(membr_potential(i));
		return out;};

private:
	vec membr_potential;			// [mV]
	vec input_train;				// [A]
	vec resistance;					// [MOhm]
	vec capacitance;				// [nF]
	mat connectivity_mat;			// [-]
	vec external;					// [mV]
	vec noise;						// [mV]
	double t_const;					// [msec]
	vec t_reset;					// [msec]
	double refract_period;			// [msec]
	const double spike_thres = -30;	// [mV]
	const double spike = 20.;		// [mV]
	const double V_rest= -70.;		// [mV]
	int N;
	ivec fired;
	uvec ind;
};
#endif /* LIF_H_ */
