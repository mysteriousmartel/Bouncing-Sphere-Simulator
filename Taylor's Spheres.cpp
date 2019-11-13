// Copyright 2019 Jennifer Campbell janeeng@bu.edu
// Copyright 2019 Haik Martirosyan haikm@bu.edu
// Copyright 2019 David Henderson dth15@bu.edu

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <math.h>

using namespace std;

// Set all floats to double, so we can do error checking easier
// Make sure to do error checking

std::vector<double> vMinus(array1, array2)
{
    finalArray[3]={};
    for(i=0;i<size(array1);i++)
    {
        finalArray[i] = array1[i]-array2[i]; 
    }
    return finalArray;
}

struct spheres () 
{
		float mass;
		float pos [3] ={};
		float vel [3] = {};
		std::string name;
		void setMass(float mass1);
		void setPos(float pos1[3]);
		void setVel(float vel1[3]);
		void setName(std::string name1);
};
void spheres::setMass(float mass1)
{
	mass=mass1;
}
void spheres::setPos(float pos1)
{
	pos = pos1;
}
void spheres::setVel(float vel1[3])
{
	vel=vel1;
}
void spheres::setName(std::string name1)
{
	name=name1;
}
float spheres::getMass()
{
	return (mass);
}
float spheres::getPos()
{
	return (pos);
}
float spheres::getVel()
{
	return (vel);
}
std::string spheres::getName()
{
	return (name);
}

struct universe()
{
	float radius;
	float max_col;
	std::vector ball_array;
	double time;
}
void universe::collideS(auto ball1,auto ball2)
{
	// have to figure out how to designate the specific
	// parts of the array
	std::vector delv = ball1vel - ball2vel
	std::vector delp = ball1pos - ball2pos
	float Radsum = ball1rad - ball2rad

	// not fully sure how this works yet, check
	// cppreference for more info

	auto A = pow(std::norm(delv),2);

	// still shaky about this vector thing too, jsyk
	// see cppreference for info on std::inner_product

	auto B = 2 * std::inner_product(std::begin(delp), std::end(delp), std::begin(delv), 0.0);

	auto c1 = std::norm(delp);
	auto c2 = pow(Radsum,2);
	auto C = c1 - c2;

	if (A == 0)
	{
		return 0;
	}
	else
	{
		auto disc = pow(B,2) - (4 * A * C);

		if (disc < 0)
		{
			return 0;
		}
		else
		{
			auto t_plus = ((-1) * B + sqrt(disc)) / (2 * A);
			auto t_minus = ((-1) * B - sqrt(disc)) / (2 * A);

			if (t_plus < 0 and t_minus < 0)
			{
				return 0;
			}
			else if (t_minus < 0)
			{
				return t_plus;
			}
			else
			{
				return t_minus;
			}
		}
	}
}

void universe::collideU(auto ball, auto radius)
{
	// need to finish this line to match ball.pos and
	// other vectors

	if (ball == {0,0,0})
	{
		return 0;
	}
	else
	{
		auto A = std::inner_product(std::begin(ballvel), std::end(ballvel), std::begin(ballvel), 0.0);
		auto B = 2 * std::inner_product(std::begin(ballpos), std::end(ballpos), std::begin(ballvel), 0.0);

		auto c1 = std::inner_product(std::begin(ballpos), std::end(ballpos), std::begin(ballpos), 0.0);
		auto c2 = pow(ballrad,2);
		auto C = c1 - c2;

		auto disc = pow(B,2) - (4 * A * C);

		if (A == 0)
		{
			return 0;
		}
		else
		{
			auto t_plus = ((-1) * B + sqrt(disc)) / (2 * A);
			auto t_minus = ((-1) * B - sqrt(disc)) / (2 * A);

			if (t_plus < 0 and t_minus < 0)
			{
				return 0;
			}
			else if (t_minus < 0)
			{
				return t_plus;
			}
			else
			{
				return t_minus;
			}
		}
	}
}

void universe::realColl(auto ball1, auto ball2, float t)
{
	std::vector<double> delv = ball1.vel - ball2.vel;
	std::vector<double> vel1_t = t * ball1.vel;
	std::vector<double> vel2_t = t * ball2.vel;
	std::vector<double> r1 = ball1.pos + vel1_t;
	std::vector<double> r2 = ball2.pos + vel2_t;
	std::vector<double> delr = r1 - r2;

	auto realCheck = std::inner_product(std::begin(delr), std::end(delr), std::begin(delv), 0.0);

	return realCheck;
}

void universe::Ucollision(auto ball1, auto ball_array, auto t, auto tot_t)
{
	auto magp=std::norm(ball1.pos);
	upos[3]= {};
	for(i=0; i<size(ball1.pos); i++)
	{
	    upos[i]=upos[i]/magp;
	}
	norm_v=std::inner_product(begin(ball1.vel), end(ball1.vel), begin(upos), 0.0);
	tan_v=minus(ball1.vel,norm_v);
	ball1.vel = minus(tan_v,norm_v);
	ball1.bounce++;
	
	std::cout<<ball1.name<<" collided with the universe\n";
	Ball::Print(ball1);
	std::cout<<"At time "<<tot_t;
	
	for(i=0;i<size(ball_array);i++)
	{
	    ball_array[i].time=t;
	}
	return 0;
}

void universe::Scollision(auto ball1, auto ball2, auto ball_array, auto t, auto tot_t)
{
	p1 = ball1.pos;
	p2 = ball2.pos;
	v1 = ball1.vel;
	v2 = ball2.vel;
	m1 = ball1.mass;
	m2 = ball2.mass;
	
	vp = vMinus(v1,v2);
	v_p = vMinus(v2,v1);
	
	rp = vMinus(p1,p2);
	r_p = vMinus(p2,p1);
	
	
	return 0;
}

void universe::update_pos(auto ball_array)
{
	return 0;
}

auto universe::energy(auto ball_array)
{
	return 0;
}

auto universe::momentum(auto ball_array)
{
	return 0;
}

auto main (auto rad, int maxCol)
{
	std::cout << "Please enter the mass, radius, x/y/z position, x/y/z velocity" << "\n";
	std::cout << "and name of each sphere"<<"\n"<<"When complete, use EOF / Ctrl-D to stop entering";
	while(eof() == 0)
	{
		
	}
}

int main()
{
    cout<<"Hello World";

    return 0;
}