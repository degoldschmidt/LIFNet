/*****************************************************************************
 *  main.cpp                                                                 *
 *                                                                           *
 *  Created on:   Oct 04, 2014                                               *
 *  Author:       Dennis Goldschmidt                                         *
 *  Email:        goldschmidtd@ini.phys.ethz.ch                              *
 *                                                                           *
 *                                                                           *
 *  Copyright (C) 2014 by Dennis Goldschmidt                                 *
 *  																		 *
 *  This file is part of the program NaviSim                                 *
 *                                                                           *
 *  NaviSim is free software: you can redistribute it and/or modify          *
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

#include <cstdio>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <armadillo>
#include "izhi.h"
#include "lif.h"
#include "timer.h"
using namespace std;
using namespace arma;

const int neurons = 100;	// number of neurons
const double T = 1000.;		// total maximum time [msec]
const double dt = 0.1;		// time interval [msec]
double t = 0.;				// time variable
int i = 0;					// discrete timesteps

Izhi * myNetwork;

ofstream data;

mat spikes;
mat potential;


double rand(double mean, double stddev){
	static random_device e{};
	static normal_distribution<double> d(mean, stddev);
	return d(e);
}

int main(){
	Timer timer(true);
	data.open("./data/results.dat");
	myNetwork = new Izhi(neurons);
	spikes.zeros(neurons, int(T/dt));
	potential.zeros(neurons, int(T/dt));

	while(t < T){
		data << t << "\t" << myNetwork->V(0) << "\t" << myNetwork->I(0)<< "\t" << myNetwork->V(1) << "\t" << myNetwork->I(1) << endl;
		if(i%1000 == 0)
			printf("t = %f secs simulated\n", t/1000.);
		//if(t>T/4.)
			//myNetwork->set_input(10.);
		myNetwork->set_noiseE(5.*rand(0., 1.));
		myNetwork->set_noiseI(2.*rand(0., 1.));
		myNetwork->update(dt);
		spikes.col(i) = myNetwork->Sp();
		potential.col(i) = myNetwork->V();
		t += dt;
		i++;
	}
	spikes.save("./data/spikes.mat", raw_ascii);
	spikes.save("./data/potential.mat", raw_ascii);

	delete myNetwork;
	data.close();
	auto elapsed_secs_cl = timer.Elapsed();
	printf("%4.3f s. Done.\n", elapsed_secs_cl.count()/1000.);
}

