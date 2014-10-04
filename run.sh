g++ main.cpp -std=c++11 -o start -O1 -larmadillo
./start
cd script/
gnuplot -persist plotx11
gnuplot -persist plotx11multi
