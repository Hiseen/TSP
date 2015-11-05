#define _AFXDLL
#ifndef _DLL
#define _DLL
#endif
#include<afxwin.h>
#include <algorithm>
#include"atlimage.h"
#include<math.h>
#include<MMSystem.h>
#include<vector>
#include<malloc.h>
#include <dshow.h>
#include <fstream> 
#include <map>
#pragma comment(lib,"strmiids.lib")
#pragma comment(lib,"quartz.lib")
using namespace std;
//需要include的
#include "hge.h"//包含hge头文件 
#include "hgeSprite.h" 
#include "hgefont.h"
#include "hgeparticle.h"
#include "hgerect.h"
#include "hgeanim.h"
#include "hgegui.h"
#include "hgeguictrls.h"
#include "menuitem.h"
#include "GfxFont.h"
#include "POINT.h"

//预编译常量

//LIB:
#pragma comment(lib,"hge.lib")
#pragma comment(lib,"hgehelp.lib")
#pragma comment(lib,"winmm.lib")
//#pragma comment(linker, "/NODEFAULTLIB:libcmt.lib")

#define MAX_Point 100
#define HEIGHT 720
#define WIDTH  1080

HGE *hge;
vector<Point> copy1;
vector<Point> box;
vector<Point> surround;
vector<Point> inter;
hgeSprite *PT;
HTEXTURE Point1;
bool start=false;
bool surring = false;
hgeFont *fot;
bool result2 = false;
bool result3 = false;
bool restart = false;
bool swi1 = false;
int type = 0;
double end1 = 0;
int endtype = 0;
struct Temp
{
	double max;
	int level;
	float a;
	float b;
	float c;
	int level2;
};
struct line
{
	bool repeat=false;
	int front;
	int back;
	double value=0;
	int ID;
	bool died;
};

struct tripack
{
	tripack();
	Point *p=NULL;
	tripack *prepack=NULL;
	Point center;
	int level;
	void Compact(tripack prepack1);
	void Update();
	bool used;
};



void tripack::Update()
{
}

tripack::tripack()
{
	p = new Point[2];
	level = -1;
}

Point GetCenter3(Point p1, Point p2, Point p3)
{
	Point p;
	p.x = (p1.x + p2.x + p3.x) / 3;
	p.y = (p1.y + p2.y + p3.y) / 3;
	p.level = -1;
	return p;
}

void tripack::Compact(tripack prepack1)
{
	if(!prepack)prepack = new tripack[2];
	for(int i = 0; i < 2;i++)
	if (prepack[i].level == -1)
	{
		prepack[i].p[0] = prepack1.p[0];
		prepack[i].p[1] = prepack1.p[1];
		prepack[i].p[2] = prepack1.p[2];
		prepack[i].level = prepack1.level;
		p[i] = prepack1.center;
	}
	else
	{
		center = GetCenter3(p[0], p[1], p[2]);
	}
}




bool ctrl = false;
//凸包内缩所要用到的函数
double direction(Point pi, Point pj, Point pk) //计算向量pkpi和向量pjpi的叉积   
{
	return (pi.x - pk.x)*(pi.y - pj.y) - (pi.y - pk.y)*(pi.x - pj.x);
}

int Distance(Point p1, Point p2)
{
	return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}

void SortPoint(vector<Point> &p, bool sort1 = true, bool sort2 = true)
{
	if (sort1)
	{
		if (sort2)sort(p.begin(), p.end(), [](const Point &p1, const Point &p2){return p1.distance > p2.distance; });
		else sort(p.begin(), p.end(), [](const Point &p1, const Point &p2){return p1.distance < p2.distance; });
	}
	p.erase(unique(p.begin(), p.end(), [](const Point &p1, const Point &p2){return p1.distance == p2.distance; }), p.end());
}

void DrawLines(vector<Point> surround,DWORD COLOR)
{

	for (int i = 0; i < surround.size(); i++)
	{
		
		fot->printf(surround[i].x, surround[i].y, HGETEXT_CENTER, "%d", i);
		if (i + 1 != surround.size())
		{
			
			hge->Gfx_RenderLine(surround[i].x, surround[i].y, surround[i + 1].x, surround[i + 1].y,COLOR);

		}
		else
			hge->Gfx_RenderLine(surround[i].x, surround[i].y, surround[0].x, surround[0].y,COLOR);
	}
}

void GetEdgePoint(vector<Point> &p, vector<Point> &np,bool reverse=true)
{
	Point p1, p2, p3, p4, p5, p6, p7, p8;
	p1 = p[0];
	p2 = p[0];
	p5 = p[0];
	p6 = p[0];
	p7 = p[0];
	p8 = p[0];
	for (int i = 0; i < p.size(); i++)
	{
		//p1 = (p1.x < p[i].x) ? p1 : p[i];
		//p2 = (p2.y < p[i].y) ? p2 : p[i];
		//p3 = (p3.x > p[i].x) ? p3 : p[i];
		//p4 = (p4.y > p[i].y) ? p4 : p[i];
		
		p5 = ((p5.y*p5.y + p5.x*p5.x) 
			< (p[i].x*p[i].x + p[i].y*p[i].y))
			? p5 : p[i];
		p6 = ((p6.y*p6.y + (WIDTH - p6.x)*(WIDTH - p6.x))
			< (p[i].y*p[i].y + (WIDTH-p[i].x)*(WIDTH-p[i].x)))
			? p6 : p[i];
		p7 = (((HEIGHT - p7.y)*(HEIGHT - p7.y) + (WIDTH - p7.x)*(WIDTH - p7.x))
			< ((HEIGHT-p[i].y)*(HEIGHT-p[i].y) + (WIDTH-p[i].x)*(WIDTH-p[i].x)))
			? p7 : p[i];
		p8 = (((HEIGHT - p8.y)*(HEIGHT - p8.y) + p8.x*p8.x) 
			< ((HEIGHT-p[i].y)*(HEIGHT-p[i].y) + p[i].x*p[i].x)) 
			? p8 : p[i];
			
	}	
	if (reverse)
	{
		SortPoint(p, true, false);
		p5 = p[0];
		p6 = p[1];
		p7 = p[2];
		p8 = p[3];
	}


	//np.push_back(p1);
	np.push_back(p5);
	//np.push_back(p2);
	np.push_back(p6);
	//np.push_back(p3);
	np.push_back(p7);
	//np.push_back(p4);
	np.push_back(p8);


	
	for (int i = 0; i < p.size(); i++)
	{
		for(int j = 0; j < np.size(); j++)
		{
			if (p[i].level == np[j].level)p[i].used = true;
		}
	}

	for (vector<Point>::iterator it = p.begin(); it != p.end();)
	{
		bool nID_in_set = it->used;
		if (nID_in_set)
			it = p.erase(it);
		else
			++it;
	}


	p.shrink_to_fit();
}

double CalculateD(vector<Point> p)
{ 
	double dis=0;
	for (int i = 0; i < p.size(); i++)
	{
		if (i + 1 != p.size())
			dis += sqrt((p[i].x - p[i + 1].x)*(p[i].x - p[i + 1].x) + (p[i].y - p[i + 1].y)*(p[i].y - p[i + 1].y));
		else
			dis += sqrt((p[i].x - p[0].x)*(p[i].x - p[0].x) + (p[i].y - p[0].y)*(p[i].y - p[0].y));
	}
	return dis;
}
/*
int AddPoint(vector<Point> &p, vector<Point> &np)
{
	
	vector<Temp> mins;
	Temp min;
	min.level = -1;
	HGE *hge = hgeCreate(HGE_VERSION);
	int pos;
    //int random = hge->Random_Int(0, p.size()-1);
	//if (p.size() <= 1)random = 0;
	
	for (int j = 0; j < p.size(); j++)
	for (int i = 0; i < np.size(); i++)
	{
		Temp min1;
		min1.level = -1;
		if (i+1 < np.size())
		{
			min1.a = (p[j].x - np[i].x)*(p[j].x - np[i].x) + (p[j].y - np[i].y)*(p[j].y - np[i].y);
			min1.b = (p[j].x - np[i + 1].x)*(p[j].x - np[i + 1].x) + (p[j].y - np[i + 1].y)*(p[j].y - np[i + 1].y);
			min1.c = (np[i].x - np[i + 1].x)*(np[i].x - np[i + 1].x) + (np[i].y - np[i + 1].y)*(np[i].y - np[i + 1].y);
			min1.level = i+1;
			min1.level2 = j;
		}
		else
		{
			
			min1.a = (p[j].x - np[i].x)*(p[j].x - np[i].x) + (p[j].y - np[i].y)*(p[j].y - np[i].y);
			min1.b = (p[j].x - np[0].x)*(p[j].x - np[0].x) + (p[j].y - np[0].y)*(p[j].y - np[0].y);
			min1.c = (np[i].x - np[0].x)*(np[i].x - np[0].x) + (np[i].y - np[0].y)*(np[i].y - np[0].y);
			min1.level2 = j;
		}
		mins.push_back(min1);
	}
	min = mins[0];
	for (int i = 0; i < mins.size(); i++)
	{
		min = ((min.a + min.b - min.c) < (mins[i].a + mins[i].b - mins[i].c)) ? min : mins[i];
	}
	if (min.level != -1)
	{
		np.insert(np.begin() + min.level, p[min.level2]);
		pos = min.level;
		np.shrink_to_fit();
	}
	else
	{
		np.insert(np.end(), p[min.level2]);
		np.shrink_to_fit();
		pos = np.size();
	}
	
	if (!p.empty())
	{
		p.erase(p.begin()+min.level2);
		p.shrink_to_fit();
	}
	return pos;

}
*/
/*
int SelectPoint(vector<Point> &p, int &a, int &b)
{
	Temp t1;
	t1.max = 0;
	for (int i = 0; i < p.size(); i++)
	{
		if (!p[i].used)
		if (i + 1 < p.size())
		{
			int a1 = Distance(p[i], p[i + 1]);
			if (a1 > t1.max)
			{
				t1.max = a1;
				t1.level = i;
				t1.level2 = i + 1;
			}
		}
		else
		{
			int a1 = Distance(p[i], p[0]);
			if (a1 > t1.max)
			{
				t1.max = a1;
				t1.level = i;
				t1.level2 = 0;
			}
		}
	}
	for (int i = 0; i < p.size(); i++)p[i].used = false;
	p[t1.level].used = true;
	p[t1.level2].used = true;
	a = t1.level;
	b = t1.level2;
	return t1.max;
}
int AddPoint(vector<Point> &p, vector<Point> &np)
{
	vector<Temp> mins;
	Temp min;
	min.level = -1;
	HGE *hge = hgeCreate(HGE_VERSION);
	int pos;
	int left, right;
	int globalc=SelectPoint(np, left, right);
	for (int i = 0; i < p.size(); i++)
	{
		Temp min1;
		min1.c = globalc;
		min1.a = (p[i].x - np[left].x)*(p[i].x - np[left].x) + (p[i].y - np[left].y)*(p[i].y - np[left].y);
		min1.b = (p[i].x - np[right].x)*(p[i].x - np[right].x) + (p[i].y - np[right].y)*(p[i].y - np[right].y);
	    min1.level = i;	
		mins.push_back(min1);
	}

	min = mins[0];
	for (int i = 0; i < mins.size(); i++)
	{
		min = ((min.a + min.b - min.c) < (mins[i].a + mins[i].b - mins[i].c)) ? min : mins[i];
	}
	    p[min.level].used = true;
		np.insert(np.begin()+right, p[min.level]);
		pos = min.level;
		np.shrink_to_fit();
		
	if (!p.empty())
	{
		p.erase(p.begin()+min.level);
		p.shrink_to_fit();
	}
	return pos;

}
*/

int AddPoint(vector<Point> &p, vector<Point> &np)
{
	vector<Temp> mins;
	Temp min;
	min.level = -1;
	HGE *hge = hgeCreate(HGE_VERSION);
	int pos;
	int random = hge->Random_Int(0, p.size() - 1);
	if (p.size() <= 1)random = 0;


	for (int i = 0; i < np.size(); i++)
	{
		Temp min1;
		min1.level = -1;
		
		if (i + 1 < np.size())
		{
			//if (((np[i].y - np[i + 1].y)*p[random].x + (np[i].x - np[i + 1].x)*p[random].y - np[i + 1].x*(np[i].y - np[i + 1].y) + np[i + 1].y*(np[i].x - np[i + 1].x))!=0)
			//{
				min1.a = (p[random].x - np[i].x)*(p[random].x - np[i].x) + (p[random].y - np[i].y)*(p[random].y - np[i].y);
				min1.b = (p[random].x - np[i + 1].x)*(p[random].x - np[i + 1].x) + (p[random].y - np[i + 1].y)*(p[random].y - np[i + 1].y);
				min1.c = (np[i].x - np[i + 1].x)*(np[i].x - np[i + 1].x) + (np[i].y - np[i + 1].y)*(np[i].y - np[i + 1].y);
				min1.level = i + 1;
			//}
		}
		else
		{
			//if (((np[i].y - np[0].y)*p[random].x + (np[i].x - np[0].x)*p[random].y - np[0].x*(np[i].y - np[0].y) + np[0].y*(np[i].x - np[0].x)) != 0)
			//{
				min1.a = (p[random].x - np[i].x)*(p[random].x - np[i].x) + (p[random].y - np[i].y)*(p[random].y - np[i].y);
				min1.b = (p[random].x - np[0].x)*(p[random].x - np[0].x) + (p[random].y - np[0].y)*(p[random].y - np[0].y);
				min1.c = (np[i].x - np[0].x)*(np[i].x - np[0].x) + (np[i].y - np[0].y)*(np[i].y - np[0].y);
			//}
		}
		mins.push_back(min1);
	}
	min = mins[0];
	for (int i = 0; i < mins.size(); i++)
	{
		min = ((min.a + min.b - min.c) < (mins[i].a + mins[i].b - mins[i].c)) ? min : mins[i];
	}
	if (min.level != -1)
	{
		np.insert(np.begin() + min.level, p[random]);
		pos = min.level;
		np.shrink_to_fit();
	}
	else
	{
		np.insert(np.end(), p[random]);
		np.shrink_to_fit();
		pos = np.size();
	}

	if (!p.empty())
	{
		p.erase(p.begin() + random);
		p.shrink_to_fit();
	}
	return pos;

}

bool ExpendPoint(vector<Point> &p, vector<Point> &np)
{
	vector<Temp> maxs;
	Temp max;
	max.max = 0;
	max.level = -1;

		for (int i = 0; i < np.size(); i++)
		{
			if (i + 1 != np.size())
			{
				for (int j = 0; j < p.size(); j++)
				{
					if ((
						((np[i + 1].y - np[i].y)*p[j].x)
						+ ((np[i].x - np[i + 1].x)*p[j].y)
						+
						(np[i + 1].x*np[i].y - np[i].x*np[i + 1].y))> 0)
					{
						Temp m1;
						m1.max = (((np[i + 1].y - np[i].y)*p[j].x)
							+ ((np[i].x - np[i + 1].x)*p[j].y)
							+ (np[i + 1].x*np[i].y - np[i].x*np[i + 1].y))
							/ (sqrt((np[i + 1].y - np[i].y)*(np[i + 1].y - np[i].y)
							+ (np[i].x - np[i + 1].x)*(np[i].x - np[i + 1].x)));

						m1.level = p[j].level;
						maxs.push_back(m1);
					}
				}
				for (int j = 0; j < maxs.size(); j++)
				{
					max = (max.max>maxs[j].max) ? max : maxs[j];
				}
				for (int j = 0; j < p.size(); j++)
				{
					if (p[j].level == max.level)
					{
						np.insert(np.begin() + i + 1, p[j]);
						p.erase(p.begin() + j);
						p.shrink_to_fit();
						np.shrink_to_fit();
						maxs.clear();
						maxs.shrink_to_fit();
						return false;
					}
				}
			}
			else
			{
				for (int j = 0; j < p.size(); j++)
				{
					if ((
						((np[0].y - np[i].y)*p[j].x)
						+ ((np[i].x - np[0].x)*p[j].y)
						+ (np[0].x*np[i].y - np[i].x*np[0].y))>0)
					{
						Temp m1;
						m1.max = (((np[0].y - np[i].y)*p[j].x)
							+ ((np[i].x - np[0].x)*p[j].y)
							+ (np[0].x*np[i].y - np[i].x*np[0].y))
							/ (sqrt((np[0].y - np[i].y)*(np[0].y - np[i].y)
							+ (np[i].x - np[0].x)*(np[i].x - np[0].x)));
						m1.level = p[j].level;
						maxs.push_back(m1);
					}
				}
				for (int j = 0; j < maxs.size(); j++)
				{
					max = (max.max>maxs[j].max) ? max : maxs[j];
				}
				for (int j = 0; j < p.size(); j++)
				{
					if (p[j].level == max.level)
					{
						np.insert(np.end(), p[j]);
						p.erase(p.begin() + j);
						p.shrink_to_fit();
						np.shrink_to_fit();
						maxs.clear();
						maxs.shrink_to_fit();
						return false;
					}
				}
			}
		}
	if (max.level == -1)return true;
}

void Expend(vector<Point> &p, vector<Point> &np)
{
	bool tor;
	do
	tor = ExpendPoint(p, np);
	while (!tor);
}


void SortAndOutFour(vector<Point> p, Point &p1, Point &p2, Point &p3, Point &p4)
{
	SortPoint(p);
	p1 = p[0];
	p2 = p[1];
	p3 = p[2];
	p4 = p[3];
}

void opt(vector<Point> &p)
{
	vector<Point> px;
	for (int i = 0; i < p.size(); i++)
	{
		/*
		if (p[i].Isnew)
		{
			for (int j = 0; j < p.size(); j++)
			{
				p[j].distance2 = Distance(p[i], p[j]);
				if (p[j].distance2 == 0)p[j].distance2 = 2147483640;
			}
		}
		SortAndOutFour(p, p1, p2, p3, p4);
		px.push_back(p1);
		px.push_back(p2);
		px.push_back(p3);
		px.push_back(p4);
		*/
		if (i + 3 < p.size())
		{
			px.push_back(p[i]);
			px.push_back(p[i + 1]);
			px.push_back(p[i + 2]);
			px.push_back(p[i + 3]);
			int dis1 = CalculateD(px);
			swap(px[1], px[2]);
			int dis2 = CalculateD(px);
			px.clear();
			if (dis2 < dis1)
				swap(p[i + 1], p[i + 2]);	
		}
		if (i + 3 == p.size())
		{
			px.push_back(p[i]);
			px.push_back(p[i + 1]);
			px.push_back(p[i + 2]);
			px.push_back(p[0]);
			int dis1 = CalculateD(px);
			swap(px[1], px[2]);
			int dis2 = CalculateD(px);
			px.clear();
			if (dis2 < dis1)
				swap(p[i + 1], p[i + 2]);
		}
	}
}



void Finalopt(vector<Point> &p)
{
	for (int i = 0; i < p.size(); i++)
		for (int j = 0; j < p.size(); j++)
		{
			if (i + 1 < p.size() && j + 1 < p.size())
			{
				double a = direction(p[i], p[i + 1], p[j]);
				double b = direction(p[i], p[i + 1], p[j + 1]);
				double c = direction(p[j], p[j + 1], p[i]);
				double d = direction(p[j], p[j + 1], p[i + 1]);
				if (a*b < 0 && c*d < 0)
				{

					int bigger = max(i, j);
					int smaller = min(i, j);
					if ((bigger - smaller) % 2)
					{
						while (bigger - smaller != 1)
						{
							swap(p[bigger], p[smaller + 1]);
							bigger--;
							smaller++;
						}
					}
					else
					{
						while (bigger != smaller)
						{
							swap(p[bigger], p[smaller + 1]);
							bigger--;
							smaller++;
						}
					}
				}
			}
			if (i +1 == p.size())
			{
				double a = direction(p[i], p[0], p[j]);
				double b = direction(p[i], p[0], p[j + 1]);
				double c = direction(p[j], p[j + 1], p[i]);
				double d = direction(p[j], p[j + 1], p[0]);
				if (a*b < 0 && c*d < 0)
				{

					int bigger = i;
					int smaller = j;
					if ((bigger - smaller) % 2)
					{
						while (bigger - smaller != 1)
						{
							swap(p[bigger], p[smaller + 1]);
							bigger--;
							smaller++;
						}
					}
					else
					{
						while (bigger != smaller)
						{
							swap(p[bigger], p[smaller + 1]);
							bigger--;
							smaller++;
						}
					}
				}

			}
		}
	

}
/*
void Myopt(vector<Point> &p)
{
	for (int i = 0; i < p.size(); i++)
	for (int j = 0; j < p.size(); j++)
	{
		if (i + 1 < p.size() && max(p[i].x, p[i + 1].x)>=p[j].x>=min(p[i].x, p[i + 1].x) && max(p[i].y, p[i + 1].y)>=p[j].y>=min(p[i].y, p[i + 1].y))
		{

			int a = (p[i].y - p[i + 1].y)*p[j].x + (p[i].x - p[i + 1].x)*p[j].y - p[i + 1].x*(p[i].y - p[i + 1].y) + p[i + 1].y*(p[i].x - p[i + 1].x);
			if (a == 0)
			{
				Point temp=p[j];
				p.erase(p.begin() + j);
				p.shrink_to_fit();
				p.insert(p.begin() + i + 1, temp);
				p.shrink_to_fit();
			}
		}
		else
		{
			double a = direction(p[j], p[i], p[0]);
			if (a == 0)
			{
				Point temp = p[j];
				p.erase(p.begin() + j);
				p.shrink_to_fit();
				p.insert(p.end(), temp);
				p.shrink_to_fit();

			}
		}
		
	}
}
*/
Point GetCenter(vector <Point> p)
{
	Point p1;
	double dis=0;
	for (int i = 0; i < p.size();i++)
	for (int j = 0; j < p.size(); j++)
	{
		double d1 = sqrt((p[i].x - p[j].x)*(p[i].x - p[j].x) + (p[i].y - p[j].y)*(p[i].y - p[j].y));
		if (dis < d1)
		{
			p1.x = (p[i].x + p[j].x) / 2;
			p1.y = (p[i].y + p[j].y) / 2;
			dis = d1;
		}
	}
	return p1;
}

void SetDistance(vector<Point> &p,Point p1)
{
	for (int i = 0; i < p.size(); i++)
		p[i].GetDis(p1.x, p1.y);
}
//从内部扩展开来用
void GetPointInside(vector<Point> &p,vector<Point> &np)
{

	SortPoint(p,true,false);
	for (int i = 0; i < 4; i++)
	{
		np.push_back(p[0]);
		p.erase(p.begin());
		p.shrink_to_fit();
	}
	Finalopt(np);
}

void GetLines(vector<line> &lines, vector<Point> p)
{
	for (int i = 0; i < p.size(); i++)
	{
		line l1;
		l1.repeat = false;
		l1.value = 0;
		l1.ID = lines.size();
		l1.died = false;
		if (i + 1 < p.size())
		{
			l1.front = p[i].level;
			l1.back = p[i + 1].level;	
		}
		else
		{
			l1.front = p[i].level;
			l1.back = p[0].level;
		}
		if(l1.front!=l1.back)
			lines.push_back(l1);
	}
	double a = 10/CalculateD(p);
	for (int i = 0; i < lines.size(); i++)
	{
		if (lines[i].value == 0)
			lines[i].value = a;
	}
	
	for (int i = 0; i < lines.size(); i++)
	for (int j = 0; j < lines.size(); j++)
	{
		if (lines[j].ID!=lines[i].ID && !lines[j].repeat && !lines[i].repeat && lines[i].front == lines[j].front && lines[j].back == lines[i].back)
		{
			lines[j].repeat = true;
			lines[i].value += lines[j].value;
		}
		if (lines[j].ID != lines[i].ID && !lines[j].repeat && !lines[i].repeat && lines[i].back == lines[j].front && lines[j].back == lines[i].front)
		{
			lines[j].repeat = true;
			lines[i].value += lines[j].value;
		}
	}
	for (vector<line>::iterator it = lines.begin(); it != lines.end();)
	{
		bool nID_in_set = it->repeat;
		if (nID_in_set)
			it = lines.erase(it);
		else
			++it;
	}
	
}


bool IsRound(vector<line> l)
{
	for (int i = 0; i < l.size(); i++)
	if (i + 1 < l.size())
	{
		if (l[i].back != l[i + 1].front)
			return false;
	}
	else
	{
		if (l[i].back != l[0].front)
			return false;
	}
	return true;
}

double FormRoute(vector<line> L, vector<Point> p)
{
	sort(L.begin(), L.end(), [](const line &L1, const line &L2){return L1.value >= L2.value; });
	vector<line> temp;

	for (int i = 0; i < L.size(); i++)
	for (int j = 0; j < L.size(); j++)
	{
		if (L[j].ID != L[i].ID && L[i].back == L[j].back)
		{
			if (L[i].value >= L[j].value)L[j].died = true;
			else L[i].died = true;
		}
		if (L[j].ID != L[i].ID && L[i].front == L[j].front)
		{
			if (L[i].value >= L[j].value)L[j].died = true;
			else L[i].died = true;
		}
	}
	for (vector<line>::iterator it = L.begin(); it != L.end();)
	{
		bool nID_in_set = it->died;
		if (nID_in_set)
		{
			
			it = L.erase(it);
		}
		else
			++it;
	}
	/*
	sort(temp.begin(), temp.end(), [](const line &L1, const line &L2){return L1.value >= L2.value; });
	

	for (int i = 0; i < L.size(); i++)
	{
		if (i + 1 < L.size() && L[i].back != L[i + 1].front)
			{
				
				for (int j = 0; j < L.size(); j++)
				{
					if (j>=i)
					{
						if (L[i].back == L[j].front)
						{
							swap(L[i + 1], L[j]);
							break;
						}
						
					}
					else
					{
						if (i-j!=1 && L[i].front == L[j].back)
						{
							for (int k = 0; k < temp.size(); k++)
							{
								if (L[j].back == temp[k].back)
								{
									swap(L[j], temp[k]);
									break;
								}
							}
						}
					}
				}
			}
	}
	if (!IsRound(L))
	{
		for (int i = 0; i < L.size(); i++)
		{
			int count = 0;
			for (int j = i; j < L.size(); j++)
			{
				if (L[i].back != L[j].front)
					count++;
			}
			if (count == L.size() - i)
			{
				line l1;
				l1.back = L[i + 1].front;
				l1.front = L[i].back;
				L.insert(L.begin() + i + 1, l1);
			}
		}
	}
	*/

	vector<Point> finalroute;
	if (!L.empty())
	{
		temp.push_back(L[0]);
		L.erase(L.begin());
	}
	while (!IsRound(temp))
	{
		for (int i = 0; i < L.size(); i++)
		{
			if (temp[temp.size() - 1].back == L[i].front)
			{
				temp.push_back(L[i]);
				L[i].died = true;
				break;
			}
			if (temp[0].front == L[i].back)
			{
				temp.insert(temp.begin(),L[i]);
				L[i].died = true;
				break;
			}
		}
		for (vector<line>::iterator it = L.begin(); it != L.end();)
		{
			bool nID_in_set = it->died;
			if (nID_in_set)
			{
				it = L.erase(it);
				break;
			}
			else
				++it;
		}
		if (temp.size() == p.size() - 1)
		{
			line l1;
			l1.back = temp[0].front;
			l1.front = temp[temp.size()-1].back;
			temp.insert(temp.end(), l1);
		}
	}
	for (int i = 0; i < temp.size(); i++)
	{
		for (int j = 0; j < p.size();j++)
		if (p[j].level == temp[i].front)
		finalroute.push_back(p[j]);
	}


	if (!IsRound(temp))
	{
		return 0;
	}
	else
	{
		return CalculateD(finalroute);
	}
	
}






void FirstTriPacking(vector<tripack> &tp, vector<Point> p)
{
	while (p.size() >= 3)
	{
		tripack *pack;
		pack = new tripack;
		pack->level = 0;
		vector<Temp> temp1;
		for (int i = 0; i < p.size(); i++)
		{
			Temp t;
			t.level = Distance(p[0], p[i]);
			t.level2 = i;
			temp1.push_back(t);
		}
		pack->p[0] = p[0];
		sort(temp1.begin(), temp1.end(), [](const Temp &p1, const Temp &p2){return p1.level < p2.level; });
		pack->p[1] = p[temp1[0].level2];
		pack->p[2] = p[temp1[1].level2];
		p[0].used = true;
		p[temp1[0].level2].used = true;
		p[temp1[1].level2].used = true;
		pack->center = GetCenter3(pack->p[0], pack->p[1], pack->p[2]);
		//删除三个被取出的点
		for (vector<Point>::iterator it = p.begin(); it != p.end();)
		{
			bool nID_in_set = it->used;
			if (nID_in_set)
				it = p.erase(it);
			else
				++it;
		}
		p.shrink_to_fit();
		tp.push_back(*pack);
	}
	int a = p.size();
	tripack *final1;
	final1 = new tripack;
	final1->level = 0;
	if (a == 1)
	{
		final1->p[0] = p[0];
	}
	if (a == 2)
	{
		final1->p[0] = p[0];
		final1->p[1] = p[1];
	}
	tp.push_back(*final1);
}

void CirclePacking(vector<tripack> &tp)
{
	while (tp.size() > 2)
	{
		vector<Temp> ints;
		for (int i = 1; i < tp.size() - 1; i++)
		{
			 Temp t1;
			 t1.level= Distance(tp[0].center, tp[i].center);
			 t1.level2 = i;
			ints.push_back(t1);
		}
		sort(ints.begin(), ints.end(), [](const Temp &p1, const Temp &p2){return p1.level < p2.level; });
		tripack *pack;
		pack = new tripack;
		for (int i = 0; i < 2; i++)
		{
			pack->Compact(tp[ints[i].level2]);
			tp[ints[i].level2].used = true;
		}
		pack->Compact(tp[0]);
		tp[0].used = true;
		for (vector<tripack>::iterator it = tp.begin(); it != tp.end();)
		{
			bool nID_in_set = it->used;
			if (nID_in_set)
				it = tp.erase(it);
			else
				++it;
		}
		tp.shrink_to_fit();
		tp.push_back(*pack);
	}

}

vector<line> lines;
int res = 0;
bool result=false;
vector<tripack> tp;
bool RenderFunc()
{




	

	hge->Gfx_BeginScene();//开始渲染 
	hge->Gfx_Clear(0xFF000000);//以某颜色清屏，OxFF000000为透明度为0的黑色
	if (start)
	{
		for (int i = 0; i < box.size(); i++)
		{
			PT->Render(box[i].x, box[i].y);
		}
		
	
	//fot->printf(100, 300, HGETEXT_CENTER, "%d", maxs.size());
	DrawLines(surround, ARGB(255,255, 0, 0));
	DrawLines(inter, ARGB(255, 0, 0, 255));
		
		fot->printf(100, 0, HGETEXT_CENTER, "%d", res);
		fot->printf(200, 0, HGETEXT_CENTER, "%d", endtype);
		fot->printf(300, 0, HGETEXT_CENTER, "%d",lines.size());
		fot->printf(400, 0, HGETEXT_CENTER, "%fl", end1);
	}
	//hge->Gfx_RenderLine(0, 0, 100, 100);
	hge->Gfx_EndScene();//结束渲染 


	
	return false;//程序正常时总是返回false，返回true将从System_Start往下执行 
} 





bool FrameFunc() 
{   
	if (!start)
	{
		box.clear();
		copy1.clear();
		int x[] = { 
			1380,
			2848,
			3510,
			457,
			3888,
			984,
			2721,
			1286,
			2716,
			738,
			1251,
			2728,
			3815,
			3683,
			1247,
			123,
			1234,
			252,
			611,
			2576,
			928,
			53,
			1807,
			274,
			2574,
			178,
			2678,
			1795,
			3384,
			3520,
			1256,
			1424,
			3913,
			3085,
			2573,
			463,
			3875,
			298,
			3479,
			2542,
			3955,
			1323,
			3447,
			2936,
			1621,
			3373,
			1393,
			3874,
			938,
			3022,
			2482,
			3854,
			376,
			2519,
			2945,
			953,
			2628,
			2097,
			890,
			2139,
			2421,
			2290,
			1115,
			2588,
			327,
			241,
			1917,
			2991,
			2573,
			19,
			3911,
			872,
			2863,
			929,
			839,
			3893,
			2178,
			3822,
			378,
			1178,
			2599,
			3416,
			2961,
			611,
			3113,
			2597,
			2586,
			161,
			1429,
			742,
			1625,
			1187,
			1787,
			22,
			3640,
			3756,
			776,
			1724,
			198,
			3950
		};
		int y[] = { 939,
			96,
			1671,
			334,
			666,
			965,
			1482,
			525,
			1432,
			1325,
			1832,
			1698,
			169,
			1533,
			1945,
			862,
			1946,
			1240,
			673,
			1676,
			1700,
			857,
			1711,
			1420,
			946,
			24,
			1825,
			962,
			1498,
			1079,
			61,
			1728,
			192,
			1528,
			1969,
			1670,
			598,
			1513,
			821,
			236,
			1743,
			280,
			1830,
			337,
			1830,
			1646,
			1368,
			1318,
			955,
			474,
			1183,
			923,
			825,
			135,
			1622,
			268,
			1479,
			981,
			1846,
			1806,
			1007,
			1810,
			1052,
			302,
			265,
			341,
			687,
			792,
			599,
			674,
			1673,
			1559,
			558,
			1766,
			620,
			102,
			1619,
			899,
			1048,
			100,
			901,
			143,
			1605,
			1384,
			885,
			1830,
			1286,
			906,
			134,
			1025,
			1651,
			706,
			1009,
			987,
			43,
			882,
			392,
			1642,
			1810,
			1558};
		
		for (int i = 0; i < MAX_Point; i++)
		{
			//x[i] /= 10;
			//y[i] /= 10;
			Point *p=new Point;
			//p->Input(x[i], y[i],i);
			p->init(i);
			box.push_back(*p);
			copy1.push_back(*p);
		}

		SortPoint(box);
		start = true;
	}

	if (start)
	{
		//凸包收缩
		if (hge->Input_KeyUp(HGEK_1))type = 1;
		if (hge->Input_KeyUp(HGEK_2))type = 2;
		if (hge->Input_KeyUp(HGEK_3))type = 3;
		if (hge->Input_KeyUp(HGEK_4))type = 4;
		if (hge->Input_KeyUp(HGEK_5))restart = true;

		if (type == 4)
		{
			if (surround.empty())
			{

				GetEdgePoint(box, surround, false);
				Expend(box, surround);
				Finalopt(surround);
				Point p1 = GetCenter(surround);
				SetDistance(box, p1);
				FirstTriPacking(tp, box);
				CirclePacking(tp);
			}

		}
		if (type == 2){
			if (surround.empty())
			{

				GetEdgePoint(box, surround, false);
				Expend(box, surround);
				Finalopt(surround);
				Point p1 = GetCenter(surround);
				SetDistance(box, p1);
				
				
			}
			if (!box.empty())

			{
				AddPoint(box, surround);
				Finalopt(surround);
			}
			else
				Expend(box, surround);
		}
		if (type == 1)
		{
			if (surround.empty())
			{

				GetEdgePoint(box, surround, false);
				Expend(box, surround);
				Finalopt(surround);
				Point p1=GetCenter(surround);
				SetDistance(box, p1);
				SortPoint(box);
			}
			if (inter.empty() && !box.empty() && box.size() >= 4)
			{
				GetEdgePoint(box, inter, false);
				Expend(box, inter);
				//SortPoint(inter,false);
			}
			else
			{
				if (!inter.empty())
				{
					AddPoint(inter, surround);
					Finalopt(surround);
				}
			}
			if (inter.empty() && !box.empty() && box.size() < 4)
			{
				    AddPoint(box, surround);
					//opt(surround);
					Finalopt(surround);
			}

		}
		if (box.empty() && inter.empty() && hge->Input_KeyUp(HGEK_ENTER))
		{
			Finalopt(surround);
		}
		if (type == 3)
		{
			if (surround.empty())
			{
				GetEdgePoint(box, inter, false);
				Expend(box, inter);
				Finalopt(inter);
				Point p1 = GetCenter(inter);
				SetDistance(box, p1);
				GetPointInside(box, surround);	
			}

			if (!box.empty())AddPoint(box, surround);
			else if(!inter.empty())AddPoint(inter, surround);
			Finalopt(surround);
		}
	}
	/*
	if (inter.empty() && box.empty())
	{
		
		restart = true;
		double res = CalculateD(surround);
		if (end1 == 0)end1 = res;
		if (res<end1)
		{
			end1 = res;
			endtype = type;
		}
		type++;
	}
	if (type > 3)type = 2;
*/
	if (restart)
	{
		//GetLines(lines, surround);
		//res=FormRoute(lines, surround);
		box.resize(copy1.size());
		memcpy(&box[0], &copy1[0], copy1.size() * sizeof(Point));
		surround.clear();
		inter.clear();
		restart = false;
		ctrl = false;
		result2 = false;
	}
	
	return false;//总是返回false 
} 


int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)//WinMain函数，程序的入口。 
{ 
hge=hgeCreate(HGE_VERSION);//使用hgeCreate函数创建HGE接口，参数必须传递正确的HGE_VERSION,它是在hge.h中定义的 
hge->System_SetState(HGE_SCREENWIDTH, 1080);//屏幕宽度 
hge->System_SetState(HGE_SCREENHEIGHT,720);//屏幕高度
hge->System_SetState(HGE_FRAMEFUNC, FrameFunc);//设置逻辑函数为FrameFunc函数 
hge->System_SetState(HGE_RENDERFUNC,RenderFunc);//设置绘制函数为RenderFunc函数 
hge->System_SetState(HGE_TITLE, "TSP");//设置窗口标题 
hge->System_SetState(HGE_WINDOWED,true);//设置使用窗口模式 
hge->System_SetState(HGE_LOGFILE,"game.log.txt");
hge->System_SetState(HGE_USESOUND,true);
hge->System_SetState(HGE_ZBUFFER,true); 
hge->System_SetState(HGE_SHOWSPLASH,false);
hge->System_SetState(HGE_SCREENBPP, 32);
hge->System_SetState(HGE_FPS,300);
hge->System_SetState(HGE_DONTSUSPEND, true);
hge->System_SetState(HGE_HIDEMOUSE,true);
//hge->System_SetState(HGE_FXVOLUME, 0);
if(hge->System_Initiate())//用hge类的System_Initiate()方法，检测初始化是否有错误出现。 
{ 

	Point1 = hge->Texture_Load("Point.png");
	PT = new hgeSprite(Point1, 0, 0, hge->Texture_GetWidth(Point1), hge->Texture_GetWidth(Point1));
	PT->SetHotSpot(PT->GetWidth() / 2, PT->GetHeight() / 2);
	PT->SetZ(0.1);
	fot = new hgeFont("font1.fnt");
	fot->SetZ(0.1);
	fot->SetScale(0.5);


	hge->System_Start();





} 

//各种释放（载入了就一定要释放！！！)

hge->System_Shutdown();
hge->Release();
return 0;
} 