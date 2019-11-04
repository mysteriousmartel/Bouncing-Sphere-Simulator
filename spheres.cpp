#include <iostream>
using namespace std;

class spheres () 
{
		float mass;
		float pos [3] ={};
		float vel [3] = {};
		std::string name;
		void setMass( float mass1);
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





int main (float rad, int maxCol)
{
	std::cout << "Please enter the mass, radius, x/y/z position, x/y/z velocity" << "\n";
	std::cout << "and name of each sphere"<<"\n"<<"When complete, use EOF / Ctrl-D to stop entering";
	while(eof() == 0)
	{
		
	}
}