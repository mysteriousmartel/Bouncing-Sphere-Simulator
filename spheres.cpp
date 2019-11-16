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

	Vec(double x1, double x2, double x3)
	{
		vec.push_back(x1);
		vec.push_back(x2);
		vec.push_back(x3);
	}

	double operator[](int i)
	{
		return vec.at(i);
	}

	Vec operator+(Vec vec2)
	{
		double x1 = vec[0] + vec2.vec[0];
		double x2 = vec[1] + vec2.vec[1];
		double x3 = vec[2] + vec2.vec[2];
		Vec vout(x1,x2,x3);
		return vout;
	}


	double Norm()
	{
		double vout = sqrt(pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2));
		return vout;
	}

	double operator*(Vec vec2)
	{
		double vout = vec[0] * vec2.vec[0] + vec[1] * vec2.vec[1] + vec[2] * vec2.vec[1];
		return vout;
	}

	Vec scalarMult(double c)
	{
		auto x1 = c * vec[0];
		auto x2 = c * vec[1];
		auto x3 = c * vec[2];
		Vec vout(x1, x2, x3);
		return vout;
	}


};

struct Ball 
{
		auto mass, rad;
		Vec pos, vel;
		string name;
		auto time, bounce;

		Ball(){}

		Ball(double mass, double rad, Vec pos, Vec vec, string name)
		{
			this->mass = mass;
			this->rad = rad;
			this->pos = pos;
			this->vel = vel;
			this->name = name;
			time = 0.0;
			bounce = 0;
		}

		void updatePos(Ball ball, auto t)
		{
			ball.pos = ball.pos + ball.vel.scalarMult(t);
		}
};

struct universe
{
	auto radius;
	auto max_col;
	auto time;

	auto collideS(Ball ball1, Ball ball2)
	{
		Vec delv = ball1.vel + ball2.vel.scalarMult(-1);
		Vec delp = ball1.pos + ball2.pos.scalarMult(-1);
		auto Radsum = ball1.rad - ball2.rad;
		auto A = pow(delv.Norm(),2);
		auto B = 2 * (delp * delv);
		auto c1 = delp.Norm();
		auto c2 = pow(Radsum,2);
		auto C = c1 - c2;

		if (A == 0)
		{
			return false;
		}
		else
		{
			auto disc = pow(B,2) - (4 * A * C);

			if (disc < 0)
			{
				return false;
			}
			else
			{
				auto t_plus = ((-1) * B + sqrt(disc)) / (2 * A);
				auto t_minus = ((-1) * B - sqrt(disc)) / (2 * A);

				if (t_plus < 0 and t_minus < 0)
				{
					return false;
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

	char collideU(Ball ball, double radius)
	{
		// unsure if this works

		if (all_of(ball.pos == {0,0,0}))
		{
			return false;
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
				return false;
			}
			else
			{
				auto t_plus = ((-1) * B + sqrt(disc)) / (2 * A);
				auto t_minus = ((-1) * B - sqrt(disc)) / (2 * A);

				if (t_plus < 0 and t_minus < 0)
				{
					return false;
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

	double realColl(Ball ball1, Ball ball2, double t)
	{
		Vec delv = ball1.vel + ball2.vel.scalarMult(-1);
		Vec vel1_t = ball1.vel.scalarMult(t);
		Vec vel2_t = ball2.vel.scalarMult(t);
		Vec r1 = ball1.pos + vel1_t;
		Vec r2 = ball2.pos + vel2_t;
		Vec delr = r1 + r2.scalarMult(-1);

		double realCheck = delr * delv;

		return realCheck;
	}

	void Ucollision(Ball ball1, vector<Ball> ball_array, double t, double tot_t)
	{
		auto magp=ball1.pos.Norm();
		Vec upos = ball1.pos.scalarMult(1/magp);

		

		norm_v= ball1.vel * upos;
		tan_v=ball1.vel + norm_v.scalarMult(-1);
		ball1.vel = tan_v + norm_v.scalarMult(-1);
		ball1.bounce++;

		cout<<"time of event: "<<tot_t<<"\n";
		
		cout<<"reflecting "<<ball1.name<<"\n";

		for(int i=0; i<ball_array.size(); i++)
	    {
	    	cout<<ball_array[i].name<<" m="<<ball_array[i].mass<<" R="<<ball_array[i].rad<<" p="
	    	<<ball_array[i].pos[0]<<","<<ball_array[i].pos[1]<<","<<ball_array[i].pos[2]<<") v=("<<
	    	ball_array[i].vel[0]<<","<<ball_array[i].vel[1]<<","<<ball_array[i].vel[2]<<") bounces="
	    	<<ball_array[i].bounce<<"\n";
	    }
		
		
		
		for(int i=0;i<size(ball_array);i++)
		{
		    ball_array[i].time=t;
		}
	}

	void Scollision(Ball ball1, Ball ball2, vector<Ball> ball_array, double t, double tot_t)
	{
		Vec p1 = ball_1.pos;
	    Vec p2 = ball_2.pos;
	    Vec v1 = ball_1.vel;
	    Vec v2 = ball_2.vel;
	    double m1 = ball_1.mass;
	    double m2 = ball_2.mass;

	    Vec del_v12 = v1+v2.scalarMult(-1);
	    Vec del_v21 =  del_p21.scalarMult(-1);
	    Vec del_p12 = p1+p2.scalarMult(-1);
	    Vec del_p21 =  del_p21.scalarMult(-1);

	    double ball1_ratio = 2 * m2/(m1+m2) * (del_v12*del_p12)/(del_p12.Norm());
	    double ball2_ratio = 2 * m1/(m1+m2) * (del_v21*del_p21)/(del_p21.Norm());

	    ball1.vel = ball1.vel + del_p12.scalarMult(-ball1_ratio);
	    ball2.vel = ball2.vel + del_p21.scalarMult(-ball2_ratio);

	    ball1.bounce++;
	    ball2.bounce++;

	    cout<<"\ntime of event: "<<tot_t<<"\n";

	    cout<<"colliding "<<ball1.name<<" "<<ball2.name<<"\n";

	    for(int i=0; i<ball_array.size(); i++)
	    {
	    	cout<<ball_array[i].name<<" m="<<ball_array[i].mass<<" R="<<ball_array[i].rad<<" p="
	    	<<ball_array[i].pos[0]<<","<<ball_array[i].pos[1]<<","<<ball_array[i].pos[2]<<") v=("<<
	    	ball_array[i].vel[0]<<","<<ball_array[i].vel[1]<<","<<ball_array[i].vel[2]<<") bounces="
	    	<<ball_array[i].bounce<<"\n";
	    }
	}

	void updatePos(vector<Ball> ball_array)
	{
		for (int i = 0; i < ball_array.size(); i++)
		{
			ball_array[i].updatePos();
			ball_array[i].time = 0.0;
		}
	}

	void energy(vector<Ball> ball_array)
	{
		tot_E=0.0;
		for(int i=0;i<ball_array.size();i++)
		{
			tot_E += (1/2)*ball_array[i].mass*pow(ball_array[i].vel.Norm(),2);
		}
		cout<<"energy: "<<tot_E<<"\n";
	}

	void momentum(vector<Ball> ball_array)
	{
		Vec mom(0.0, 0.0, 0.0);
		for(int i=0; i<ball_array.size();i++)
		{
			mom = mom + ball_array[i].vel.scalarMult(ball_array[i].mass);
		}
		cout<<"momentum: ("<<mom[0]<<", "<<mom[1]<<", "<<mom[2]<<")\n";

	}
};

int main (int argc, char **argv)
{
	double univ_rad = atof(*(argv+1));
	int univ_coll = atof(*(argv+2));

	vector<Ball> ball_array;
	double total_time = 0.0;

	double mass, radius, x0, y0, z0, vx, vy, vz;
	string name;

	cout << "Please enter the mass, radius, x/y/z position, x/y/z velocity" << "\n";
	cout << "and name of each sphere"<<"\n"<<"When complete, use EOF / Ctrl-D to stop entering";
	while(!cin.eof())
	{
		cin>>mass>>radius>>x0>>y0>>x0>>vx>>vy>>vz>>name;
		Vec pos(x0,y0,z0);
		Vec vel(vx,vy,vz);

		ball_array.push_back(Ball(mass,radius,name,pos,vel));
	}

	universe u = universe(univ_rad,univ_coll);

	auto minty = u.collideU(ball_array[0],univ_rad);
	vector<int> colliders{0,0};

	cout<<"Here are the initial conditions.\nuniverse radius "<<univ_rad<<"\nmax collisions "
	<<univ_coll<<endl;

	for(int i=0; i<ball_array.size(); i++)
	{
		cout<<ball_array[i].name<<" m="<<ball_array[i].mass<<" R="<<ball_array[i].rad<<" p="
	    <<ball_array[i].pos[0]<<","<<ball_array[i].pos[1]<<","<<ball_array[i].pos[2]<<") v=("<<
	    ball_array[i].vel[0]<<","<<ball_array[i].vel[1]<<","<<ball_array[i].vel[2]<<") bounces="
	    <<ball_array[i].bounce<<"\n";
	}
	u.energy(ball_array);
	u.momentum(ball_array);
	cout<<"\nHere are the events"<<endl;


	while (!ball_array.empty())
	{
		if (ball_array.size()==1)
		{
			minty = u.collideU(ball_array[0],univ_rad);
			colliders = {0,-1};
		}
		else
		{
			for (int i=0; i<ball_array.size()-1;i++)
			{
				double t = u.collideU(ball_array[i],univ_rad);
				if (t != false and t<minty)
				{
					if (colliders != (i,-1))
					{
						minty = t;
						colliders = (i,-1);
					}
				}
				for (int j=i+1;j<ball_array.size();j++)
				{
					t = u.collideS(ball_array[i],ball_array[j]);
					if (t != false and t < minty)
					{
						if (colliders != (i,j))
						{
							collCheck = u.realColl(ball_array[i],ball_array[j],t);
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
			u.energy(ball_array);
			u.momentum(ball_array);
			u.Ucollision(ball_array[colliders[0]],ball_array,minty,total_time);
			
			ball_array[colliders[0]].bounce = ball_array[colliders[0]].bounce + 1;

			if (ball_array[colliders[0]].bounce == univ_coll)
			{
				cout<<"disappear "<<ball_array[colliders[0]].name<<endl;
				ball_array.erase(colliders[0]);
			}
		}
		else
		{
			u.energy(ball_array);
		    u.momentum(ball_array);
			u.Scollision(ball_array[colliders[0]], ball_array[colliders[1]], ball_array, minty, totalTime);
		    ball_array[colliders[0]].bounce = ball_array[colliders[0]].bounce + 1;
		    ball_array[colliders[1]].bounce = ball_array[colliders[1]].bounce + 1;
	      if (ball_array[colliders[0]].bounce > = univ_coll and colliders[1] >= univ_coll)
	      {
	      	cout<<"disappear "<<ball_array[colliders[0]].name<<endl;
	      	cout<<"disappear "<<ball_array[colliders[1]].name<<endl;
	      	ball_array.erase(colliders[1]);
	      	ball_array.erase(colliders[0]);
	      }
	      else if (ball_array[colliders[1]].bounce >= univ_coll)
	      {
	      	cout<<"disappear "<<ball_array[colliders[1]].name<<endl;
	      	ball_array.erase(colliders[1]);
	      }
	      else if (ball_array[colliders[0]].bounce >= univ_coll)
	      {
	      	cout<<"disappear "<<ball_array[colliders[0]].name<<endl;
	      	ball_array.erase(colliders[0]);
	      }
		}
		u.updatePos(ball_array);
	}
	cout<<"Total time for all spheres to vanish: "<<total_time<<endl;
	return 0.0;
}