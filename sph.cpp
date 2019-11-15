// Copyright 2019 Jennifer Campbell janeeng@bu.edu
// Copyright 2019 Haik Martirosyan haikm@bu.edu
// Copyright 2019 David Henderson dth15@bu.edu

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <string>
#include <numeric>
#include <math.h>

using namespace std;

// Set all floats to double, so we can do error checking easier
// Make sure to do error checking

class spheres {
    public:
    double mass;
    double radius;
    string name;
    double pos[3] 
    double vel[3] 
    double time;

    Ball(double mass, double radius, string name, double *pos, double *vel){
        this->mass=mass;
        this->radius=radius;
        this->name=name;
        this->pos[0]=pos[0];
        this->pos[1]=pos[1];
        this->pos[2]=pos[2];
        this->vel[0]=vel[0];
        this->vel[0]=vel[1];
        this->vel[0]=vel[2];
        this->time=0.0;

    void updatePos(double t) {
        double dt = t - this->time;
        this->pos[0] += dt * this->vel[0]; 
        this->pos[1] += dt * this->vel[1]; 
        this->pos[2] += dt * this->vel[3]; 
        this->time = t;
    }
};

class universe {
	public:
    float radius;
	float max_col;
	vector<spheres> ball_array;
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
	std::vector<double> delv = ball1vel - ball2vel;
	std::vector<double> vel1_t = t * ball1vel;
	std::vector<double> vel2_t = t * ball2vel;
	std::vector<double> r1 = ball1pos + vel1_t;
	std::vector<double> r2 = ball2pos + vel2_t;
	std::vector<double> delr = r1 - r2;

	auto realCheck = std::inner_product(std::begin(delr), std::end(delr), std::begin(delv), 0.0);

	return realCheck;
}

void universe::Ucollision(auto ball1, auto ball_array, auto t, auto tot_t)
{
	return 0;
}

void universe::Scollision(auto ball1, auto ball2, auto ball_array, auto t, auto tot_t)
{
    vec p1 = ball_1.pos;
    vec p2 = ball_2.pos;
    vec v1 = ball_1.vel;
    vec v2 = ball_2.vel;
    double m1 = ball_1.mass;
    double m2 = ball_2.mass;

    del_v12 = v1-v2;
    del_v21 = scalarMul(-1, del_p21);
    del_p12 = p1-p2;
    del_p21 = scalarMul(-1, del_p21);

    ball1_ratio = 2 * m2/(m1+m2) * (del_v12*del_p12)/(norm(del_p12));
    ball2_ratio = 2 * m1/(m1+m2) * (del_v21*del_p21)/(norm(del_p21));

    ball1.vec = ball1.vel - scalarMul(ball1_rat,del_p12);
    ball2.vec = ball2.vel - scalarMul(ball2_rat,del_p21);

    ball1.bounce ++;
    ball2.bounce ++;
	
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
