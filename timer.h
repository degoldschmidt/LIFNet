/*****************************************************************************
 *  timer.h                                                                  *
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

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <ostream>

class Timer {
	typedef std::chrono::high_resolution_clock high_resolution_clock;
	typedef std::chrono::milliseconds milliseconds;
public:
	explicit Timer(bool run = false)
	{
		if (run)
			Reset();
	}
	void Reset()
	{
		_start = high_resolution_clock::now();
	}
	milliseconds Elapsed() const
	{
		return std::chrono::duration_cast<milliseconds>(high_resolution_clock::now() - _start);
	}
	template <typename T, typename Traits>
	friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer)
	{
		return out << timer.Elapsed().count();
	}
private:
	high_resolution_clock::time_point _start;
};



#endif /* TIMER_H_ */
