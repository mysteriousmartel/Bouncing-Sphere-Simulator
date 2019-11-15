// Copyright 2019 Jennifer Campbell janeeng@bu.edu
// Copyright 2019 Haik Martirosyan haikm@bu.edu
// Copyright 2019 David Henderson dth15@bu.edu

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <math.h>
#include <list>

using namespace std;

// Make sure to do error checking

struct Vec
{
	vector<double> vec;

	Vec(){}

	Vec(auto x1, auto x2, auto x3)
	{
		vec.push_back(x1);
		vec.push_back(x2);
		vec.push_back(x3);
	}

	Vec operator+(vec1, vec2);
	{
		auto x1 = vec1[0] + vec2[0];
		auto x2 = vec1[1] + vec2[1];
		auto x3 = vec1[2] + vec2[2];
		Vec vout(x1,x2,x3);
		return vout;
	}

	Vec operator*(Vec vec1, Vec vec2)
	{
		auto vout = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[1];
		return vout;
	}

	Vec scalarMult(auto c, Vec vec)
	{
		auto x1 = c * vec[0];
		auto x2 = c * vec[1];
		auto x3 = c * vec[2];
		Vec vout(x1, x2, x3);
		return vout;
	}

	Vec Norm(Vec vec)
	{
		auto vout = sqrt(vec * vec);
		return vout;
	}
}

struct ball() 
{
		auto mass, rad;
		Vec pos, vel;
		string name;
		auto time, bounce;

		Ball(){}

		Ball(auto mass, auto rad, Vec pos, Vec vec, string name)
		{
			this->mass = mass;
			this->rad = rad;
			this->pos = pos;
			this->vel = vel;
			this->name = name;
			time = 0.0;
			bounce = 0;
		}

		auto updatePos(Ball ball, auto t)
		{
			ball.pos = ball.pos + scalarMult(t,ball.vel);
		}
};

struct universe()
{
	auto radius;
	auto max_col;
	vector<Ball> ball_array;
	auto time;

	auto collideS(Ball ball1, Ball ball2)
	{
		Vec delv = ball1.vel - ball2.vel
		Vec delp = ball1.pos - ball2.pos
		auto Radsum = ball1.rad - ball2.rad
		auto A = pow(Norm(delv),2);
		auto B = 2 * (delp * delv);
		auto c1 = Norm(delp);
		auto c2 = pow(Radsum,2);
		auto C = c1 - c2;

		if (A == 0)
		{
			return False;
		}
		else
		{
			auto disc = pow(B,2) - (4 * A * C);

			if (disc < 0)
			{
				return False;
			}
			else
			{
				auto t_plus = ((-1) * B + sqrt(disc)) / (2 * A);
				auto t_minus = ((-1) * B - sqrt(disc)) / (2 * A);

				if (t_plus < 0 and t_minus < 0)
				{
					return False;
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

	auto collideU(Ball ball, auto radius)
	{
		// unsure if this works

		if (all_of(ball.pos == {0,0,0}))
		{
			return False;
		}
		else
		{
			auto A = ball.vel * ball.vel;
			auto B = 2 * ball.pos * ball.vel;
			auto c1 = ball.pos * ball.pos;
			auto c2 = pow(ball.rad,2);
			auto C = c1 - c2;
			auto disc = pow(B,2) - (4 * A * C);

			if (A == 0)
			{
				return False;
			}
			else
			{
				auto t_plus = ((-1) * B + sqrt(disc)) / (2 * A);
				auto t_minus = ((-1) * B - sqrt(disc)) / (2 * A);

				if (t_plus < 0 and t_minus < 0)
				{
					return False;
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

	void realColl(Ball ball1, Ball ball2, auto t)
	{
		Vec delv = ball1.vel + scalarMult(-1,ball2.vel);
		Vec vel1_t = scalarMult(t,ball1.vel);
		Vec vel2_t = scalarMult(t,ball2.vel);
		Vec r1 = ball1.pos + vel1_t;
		Vec r2 = ball2.pos + vel2_t;
		Vec delr = r1 + scalarMult(-1,r2);

		auto realCheck = delr * delv;

		return realCheck;
	}

	void Ucollision(Ball ball1, vector<Ball> ball_array, auto t, auto tot_t)
	{
		auto magp=Norm(ball1.pos);
		Vec upos[3]= {};

		for(i=0; i<size(ball1.pos); i++)
		{
		    upos[i]=upos[i]/magp;
		}

		norm_v= ball1.vel * upos;
		tan_v=ball1.vel + scalarMult(-1,norm_v);
		ball1.vel = tan_v + scalarMult(-1,norm_v);
		ball1.bounce++;
		
		cout<<ball1.name<<" collided with the universe\n";
		Ball::Print(ball1);
		cout<<"At time "<<tot_t;
		
		for(i=0;i<size(ball_array);i++)
		{
		    ball_array[i].time=t;
		}
		return 0;
	}

	void Scollision(Ball ball1, Ball ball2, vector<Ball> ball_array, auto t, auto tot_t)
	{
		p1 = ball1.pos;
		p2 = ball2.pos;
		v1 = ball1.vel;
		v2 = ball2.vel;
		m1 = ball1.mass;
		m2 = ball2.mass;
		
		vp = v1 + scalarMult(-1,v2);
		v_p = v2 + scalarMult(-1,v1);
		
		rp = p1 + scalarMult(-1,p2);
		r_p = p2 + scalarMult(-1,p1);
		
		
		return 0;
	}

	void update_pos(vector<Ball> ball_array)
	{
		for (i = 0; i < ball_array.size(); i++)
		{
			ball_array[i].updatePos();
			ball_array[i].time = 0.0;
		}
	}

	auto energy(vector<Ball> ball_array)
	{
		return 0;
	}

	Vec momentum(vector<Ball> ball_array)
	{
		return 0;
	}
};

int main (auto rad, int maxCol)
{
	auto univ_rad = rad;
	int univ_coll = maxCol;

	vector<double> initials;
	auto total_time = 0.0;

	auto mass, radius, x0, y0, z0, vx, vy, vz;
	string name;

	cout << "Please enter the mass, radius, x/y/z position, x/y/z velocity" << "\n";
	cout << "and name of each sphere"<<"\n"<<"When complete, use EOF / Ctrl-D to stop entering";
	while(!cin.eof())
	{
		cin>>mass>>radius>>x0>>y0>>x0>>vx>>vy>>vz>>name;
		Vec pos(x0,y0,z0);
		Vec vel(vx,vy,vz);

		initials.push_back(Ball(mass,radius,name,pos,vel));
	}

	void universe = Universe(univ_rad,univ_coll,initials);
	vector<Ball> ball_array = universe.ball_array;

	auto minty = universe.collideU(ball_array[0],univ_rad);
	array<int> colliders(0,0);

	while (!ball_array.empty())
	{
		if (ball_array.size()==1)
		{
			minty = universe.collideU(ball_array[0],univ_rad);
			colliders = {0,-1};
		}
		else
		{
			for (i=0; i<ball_array.size()-1;i++)
			{
				t = universe.collideU(ball_array[i],univ_rad);
				if (t != False and t<minty)
				{
					if (colliders != (i,-1))
					{
						minty = t;
						colliders = (i,-1);
					}
				}
				for (j=i+1;j<ball_array.size();j++)
				{
					t = universe.collideS(ball_array[i],ball_array[j]);
					if (t != False and t < minty)
					{
						if (colliders != (i,j))
						{
							collCheck = universe.realColl(ball_array[i],ball_array[j],t);
							if (collCheck < 0)
							{
								minty = t;
								colliders = (i,j);
							}
						}
					}
				}
			}
		}
		total_time = total_time + minty;

		if (colliders[1] == -1)
		{
			universe.Ucollision(ball_array[colliders[0]],ball_array,minty,total_time);
			universe.energy(ball_array);
			universe.momentum(ball_array);
			ball_array[colliders[0]].bounce = ball_array[colliders[0]].bounce + 1;

			if (ball_array[colliders[0]].bounce == univ_coll)
			{
				cout<<ball_array[colliders[0]].name<<" has left"<<endl;
				ball_array.erase(colliders[0]);
			}
		}
		else
		{
			universe.Scollision(ball_array[colliders[0]], ball_array[colliders[1]], ball_array, minty, totalTime);
      universe.energy(ball_array);
      universe.momentum(ball_array);
      ball_array[colliders[0]].bounce = ball_array[colliders[0]].bounce + 1
      ball_array[colliders[1]].bounce = ball_array[colliders[1]].bounce + 1

      if (ball_array[colliders[0]].bounce > = univ_coll and colliders[1] >= univ_coll)
      {
      	cout<<ball_array[colliders[0]].name<<" has left"<<endl;
      	cout<<ball_array[colliders[1]].name<<" has left"<<endl;
      	ball_array.erase(colliders[1]);
      	ball_array.erase(colliders[0]);
      }
      else if (ball_array[colliders[1]].bounce >= univ_coll)
      {
      	cout<<ball_array[colliders[1]].name<<" has left"<<endl;
      	ball_array.erase(colliders[1]);
      }
      else if (ball_array[colliders[0]].bounce >= univ_coll)
      {
      	cout<<ball_array[colliders[0]].name<<" has left"<<endl;
      	ball_array.erase(colliders[0]);
      }
		}
		universe.updatePos(ball_array);
	}
	cout<<"Total time for all spheres to vanish: "<<total_time<<endl;
}