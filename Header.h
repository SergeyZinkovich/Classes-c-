#pragma once
#include <vector>

using std::vector;

class Point;
class Base_class;
class Line;
class PolyLine;
class Circle;

class Point {
public:
	double x, y;
	Point(double x0, double y0);
};

class Base_class {
public:
	virtual double Length() = 0;
	virtual std::vector<Point> Intersect(Line) = 0;
	virtual std::vector<Point> Intersect(PolyLine) = 0;
	virtual std::vector<Point> Intersect(Circle) = 0;
	virtual ~Base_class() {};
};

class Line : public Base_class {
public:
	vector<Point> points;
	Line(vector<Point> a);
	double Length() override;
	virtual vector<Point> Intersect(Line other) override;
	virtual vector<Point> Intersect(PolyLine other) override;
	virtual vector<Point> Intersect(Circle other) override;
	virtual ~Line() {};
};

class PolyLine : public Line {
public:
	PolyLine(vector<Point> a);
	double Length();
	virtual Line Get_line(int i);
	virtual vector<Point> Intersect(Line other);
	virtual vector<Point> Intersect(PolyLine other);
	virtual vector<Point> Intersect(Circle other);
	virtual ~PolyLine() {};
};

class Circle : public Base_class {
public:
	double x, y, r;
	Circle(double x0, double y0, double r0);
	double Length();
	virtual vector<Point> Intersect(Line other);
	virtual vector<Point> Intersect(PolyLine other);
	virtual vector<Point> Intersect(Circle other);
	virtual ~Circle() {};
};

