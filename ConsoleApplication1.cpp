#define CATCH_CONFIG_MAIN
#define _USE_MATH_DEFINES
#define EPS 0.000000001
#include "Header.h"
#include <iostream>
#include <vector>
#include <math.h>
#include "catch.hpp"

using std::vector;

Point::Point(double x0, double y0) {
	x = x0;
	y = y0;
}

Line::Line(vector<Point> a) {
	points = a;
}

double Line::Length() {
	return sqrt(pow(points[1].x - points[0].x, 2) + pow(points[1].y + points[0].y, 2));
}

vector<Point> Line::Intersect(Line other) {
	vector<Point> result;
	double x1, y1, x2, y2, x3, y3, x4, y4;
	x1 = this->points[0].x;
	y1 = this->points[0].y;
	x2 = this->points[1].x;
	y2 = this->points[1].y;
	x3 = other.points[0].x;
	y3 = other.points[0].y;
	x4 = other.points[1].x;
	y4 = other.points[1].y;

	double d = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
	double u_a = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / d;
	double u_b = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / d;

	if (u_a >= 0 && u_a <= 1 && u_b >= 0 && u_b <= 1) {
		double ans_x = x1 + u_a * (x2 - x1);
		double ans_y = y1 + u_a * (y2 - y1);
		result.push_back(Point(ans_x, ans_y));
	}
	return result;
}

vector<Point> Line::Intersect(PolyLine other) {
	vector<Point> result;
	for (int i = 1; i < other.points.size(); i++) {
		Line l = other.Get_line(i);
		vector<Point> r = l.Intersect(*this);
		for (int j = 0; j < r.size(); j++) {
			result.push_back(r[j]);
		}
	}
	return result;
}

vector<Point> Line::Intersect(Circle other) {
	vector<Point> result;
	Point a = this->points[0];
	Point b = this->points[1];
	Point o = Point(other.x, other.y);
	double r = other.r;

	a.x -= o.x; a.y -= o.y;
	b.x -= o.x; b.y -= o.y;

	double A = a.y - b.y;
	double B = b.x - a.x;
	double C = a.x * b.y - b.x * a.y;

	Point q = Point(0, 0);
	q.x = -(A * C) / (A * A + B * B);
	q.y = -(B * C) / (A * A + B * B);
	if (C * C <= r * r * (A * A + B * B) + EPS) {
		if (abs(C * C - r * r * (A * A + B * B)) < EPS) {
			q.x += o.x;
			q.y += o.y;
			result.push_back(q);
		}
		else {
			double d = r * r - C * C / (A * A + B * B);
			double mult = sqrt(d / (A * A + B * B));
			result.push_back(Point(q.x + B * mult + o.x, q.y - A * mult + o.y));
			result.push_back(Point(q.x - B * mult + o.x, q.y + B * mult + o.y));
		}
	}
	return result;
}

PolyLine::PolyLine(vector<Point> a) :Line(a) {}

double PolyLine::Length() {
	double l = 0;
	for (int i = 1; i < points.size(); i++) {
		l += sqrt(pow(points[i].x - points[i - 1].x, 2) + pow(points[i].y + points[i - 1].y, 2));
	}
	return l;
}

Line PolyLine::Get_line(int i) {
	vector<Point> a;
	a.push_back(Point(points[i - 1].x, points[i - 1].y));
	a.push_back(Point(points[i].x, points[i].y));
	return Line(a);
}

vector<Point> PolyLine::Intersect(Line other) {
	return other.Intersect(*this);
}

vector<Point> PolyLine::Intersect(PolyLine other) {
	vector<Point> result;
	for (int i = 1; i < this->points.size(); i++) {
		vector<Point> r = this->Get_line(i).Intersect(other);
		for (int j = 0; j < r.size(); j++) {
			result.push_back(r[j]);
		}
	}
	return result;
}

vector<Point> PolyLine::Intersect(Circle other) {
	vector<Point> result;
	for (int i = 1; i < this->points.size(); i++) {
		vector<Point> r = this->Get_line(i).Intersect(other);
		for (int j = 0; j < r.size(); j++) {
			result.push_back(r[j]);
		}
	}
	return result;
}

Circle::Circle(double x0, double y0, double r0) {
	x = x0;
	y = y0;
	r = r0;
}

double Circle::Length() {
	return 2 * M_PI * r;
}

vector<Point> Circle::Intersect(Line other) {
	return other.Intersect(*this);
}

vector<Point> Circle::Intersect(PolyLine other) {
	return other.Intersect(*this);
}

vector<Point> Circle::Intersect(Circle other) {
	vector<Point> result;

	double d = sqrt(pow(other.x - this->x, 2) + pow(other.y - this->y, 2));
	bool nesting = abs(other.r - this->r) > d;
	bool no_intersection = d > other.r + this->r;
	if (!(nesting || no_intersection)) {
		double b = (pow(this->r, 2) - pow(other.r, 2) + pow(d, 2)) / (2 * d);
		double a = d - b;
		Point p0 = Point(0, 0);
		p0.x = other.x + a / d * (this->x - other.x);
		p0.y = other.y + a / d * (this->y - other.y);
		if (d == other.r + this->r)
			result.push_back(p0);
		else {
			double h = sqrt(pow(other.r, 2) - pow(a, 2));
			Point p3 = Point(0, 0), p4 = Point(0, 0);
			p3.x = p0.x + (this->y - other.y) * h / d;
			p3.y = p0.y - (this->x - other.x) * h / d;
			p4.x = p0.x - (this->y - other.y) * h / d;
			p4.y = p0.y + (this->x - other.x) * h / d;
			result.push_back(p3);
			result.push_back(p4);
		}
	}
	return result;
}

TEST_CASE("Intersect line with line", "[]") {
	SECTION("Simple intersection") {
		vector<Point> e, f;
		e.push_back(Point(0.0, 0.0));
		e.push_back(Point(4.0, 4.0));
		f.push_back(Point(2.0, 0.0));
		f.push_back(Point(2.0, 4.0));

		Line line1 = Line(e);
		Line line2 = Line(f);

		Point expected_point(2.0, 2.0);
		vector<Point> result = line1.Intersect(line2);

		REQUIRE(result[0].x == expected_point.x);
		REQUIRE(result[0].y == expected_point.y);
	}

	SECTION("No intersection") {
		Point a = Point(0.0, 0.0), b = Point(4.0, 4.0), c = Point(2.0, 3.0), d = Point(2.0, 4.0);
		vector<Point> e, f;
		e.push_back(Point(0.0, 0.0));
		e.push_back(Point(4.0, 4.0));
		f.push_back(c);
		f.push_back(d);

		Line line1 = Line(e);
		Line line2 = Line(f);

		vector<Point> result = line1.Intersect(line2);

		REQUIRE(result.empty() == true);
	}

	SECTION("Same point") {
		Point a = Point(0.0, 0.0), b = Point(4.0, 4.0), c = Point(0.0, 0.0), d = Point(0.0, 2.0);
		vector<Point> e, f;
		e.push_back(a);
		e.push_back(b);
		f.push_back(c);
		f.push_back(d);

		Line line1 = Line(e);
		Line line2 = Line(f);

		Point expected_point(Point(0.0, 0.0));
		vector<Point> result = line1.Intersect(line2);

		REQUIRE(result[0].x == expected_point.x);
		REQUIRE(result[0].y == expected_point.y);
	}

	SECTION("Same lines") {
		Point a = Point(0.0, 0.0), b = Point(4.0, 4.0), c = Point(0.0, 0.0), d = Point(4.0, 4.0);
		vector<Point> e, f;
		e.push_back(a);
		e.push_back(b);
		f.push_back(c);
		f.push_back(d);

		Line line1 = Line(e);
		Line line2 = Line(f);

		vector<Point> result = line1.Intersect(line2);

		REQUIRE(result.empty() == true);
	}

	SECTION("Very close, no intersection") {
		Point a = Point(0.0, 0.0), b = Point(4.0, 4.0), c = Point(2.0, 0.0), d = Point(2.0, 1.9999999999999);
		vector<Point> e, f;
		e.push_back(a);
		e.push_back(b);
		f.push_back(c);
		f.push_back(d);

		Line line1 = Line(e);
		Line line2 = Line(f);

		vector<Point> result = line1.Intersect(line2);

		REQUIRE(result.empty() == true);
	}

}

TEST_CASE("Intersect line with polyline", "[]") {
	vector<Point> points;
	points.push_back(Point(0.0, 0.0));
	points.push_back(Point(4.0, 4.0));
	points.push_back(Point(2.0, 0.0));
	points.push_back(Point(2.0, 4.0));
	PolyLine broken_line(points);

	vector<Point> a;
	a.push_back(Point(0.0, 3.0));
	a.push_back(Point(4.0, 3.0));
	Line line(a);

	vector<Point> expected_points;
	expected_points.push_back(Point(3.0, 3.0));
	expected_points.push_back(Point(3.5, 3.0));
	expected_points.push_back(Point(2.0, 3.0));
	vector<Point> result = line.Intersect(broken_line);
	REQUIRE(result.size() == expected_points.size());
	for (int i = 0; i < result.size(); i++) {
		REQUIRE(result[i].x == expected_points[i].x);
		REQUIRE(result[i].y == expected_points[i].y);
	}
}

TEST_CASE("Intersect polyline with polyline") {
	vector<Point> points;
	points.push_back(Point(6.0, 4.0));
	points.push_back(Point(0.0, 1.0));
	points.push_back(Point(5.0, 1.0));
	PolyLine first_bline(points);

	points.clear();
	points.push_back(Point(1.0, 3.0));
	points.push_back(Point(4.0, 0.0));
	points.push_back(Point(4.0, 4.0));
	PolyLine second_bline(points);

	vector<Point> expected_points;
	expected_points.push_back(Point(2.0, 2.0));
	expected_points.push_back(Point(4.0, 3.0));
	expected_points.push_back(Point(3.0, 1.0));
	expected_points.push_back(Point(4.0, 1.0));

	vector<Point> result = first_bline.Intersect(second_bline);
	REQUIRE(result.size() == expected_points.size());
	for (int i = 0; i < result.size(); i++) {
		REQUIRE(result[i].x == expected_points[i].x);
		REQUIRE(result[i].y == expected_points[i].y);
	}
}

TEST_CASE("Intersect line with circle") {
	SECTION("In two points") {
		Circle circle(2.0, 2.0, 2);

		vector<Point> b;
		b.push_back(Point(1.0, 5.0));
		b.push_back(Point(5.0, 1.0));
		Line line(b);

		vector<Point> expected_points;
		expected_points.push_back(Point(4.0, 2.0));
		expected_points.push_back(Point(2.0, 4.0));

		vector<Point> result = line.Intersect(circle);
		REQUIRE(result.size() == 2);
		for (int i = 0; i < result.size(); i++) {
			REQUIRE(expected_points[i].x == result[i].x);
			REQUIRE(expected_points[i].y == result[i].y);
		}
	}

	SECTION("No intersection") {
		Circle circle(2.0, 2.0, 2);

		vector<Point> b;
		b.push_back(Point(5.0, 5.0));
		b.push_back(Point(5.0, 1.0));
		Line line(b);

		vector<Point> result = circle.Intersect(line);
		REQUIRE(result.size() == 0);
	}

	SECTION("In one point") {
		Circle circle(2.0, 2.0, 2);

		vector<Point> b;
		b.push_back(Point(4.0, 5.0));
		b.push_back(Point(4.0, 1.0));
		Line line(b);

		Point expected_point(Point(4.0, 2.0));

		vector<Point> result = line.Intersect(circle);
		REQUIRE(result.size() == 1);
		REQUIRE(expected_point.x == result[0].x);
		REQUIRE(expected_point.y == result[0].y);
	}
}

TEST_CASE("Intersect polyline line with circle") {
	vector<Point> points;
	points.push_back(Point(2.0, 5.0));
	points.push_back(Point(2.0, -1.0));
	points.push_back(Point(4.0, -1.0));
	points.push_back(Point(4.0, 5.0));
	PolyLine broken_line(points);

	Circle circle(2.0, 2.0, 2.0);

	vector<Point> expected_points;
	expected_points.push_back(Point(2.0, 0.0));
	expected_points.push_back(Point(2.0, 2.0));
	expected_points.push_back(Point(4.0, 2.0));

	vector<Point> result = circle.Intersect(broken_line);
	REQUIRE(result.size() == 3);
	for (int i = 0; i < result.size(); i++) {
		REQUIRE(expected_points[i].x == result[i].x);
		REQUIRE(expected_points[i].y == result[i].y);
	}
}

TEST_CASE("Intersect circle with circle", "[]") {
	SECTION("Intersect in two points") {
		Circle first_circle(0.0, 0.0, 2.0);

		Circle second_circle(2.0, 2.0, 2.0);

		vector<Point> expected_points;
		expected_points.push_back(Point(0.0, 2.0));
		expected_points.push_back(Point(2.0, 0.0));

		vector<Point> result = first_circle.Intersect(second_circle);
		REQUIRE(result.size() == expected_points.size());
		for (int i = 0; i < result.size(); i++) {
			REQUIRE(result[i].x == expected_points[i].x);
			REQUIRE(result[i].y == expected_points[i].y);
		}
	}

	SECTION("Intersect in one point") {
		Circle first_circle(0.0, 0.0, 3.0);

		Circle second_circle(0.0, 6.0, 3.0);

		vector<Point> expected_points;
		expected_points.push_back(Point(0.0, 3.0));

		vector<Point> result = first_circle.Intersect(second_circle);
		REQUIRE(result.size() == expected_points.size());
		for (int i = 0; i < result.size(); i++) {
			REQUIRE(result[i].x == expected_points[i].x);
			REQUIRE(result[i].y == expected_points[i].y);
		}
	}

	SECTION("No intersection") {
		Circle first_circle(0.0, 0.0, 3.0);

		Circle second_circle(5.0, 0.0, 1.0);

		vector<Point> result = first_circle.Intersect(second_circle);
		REQUIRE(result.size() == 0);
		int kek;
	}

	SECTION("Nested circles") {
		Circle first_circle(3.0, 3.0, 3.0);

		Circle second_circle(3.0, 3.0, 1.0);

		vector<Point> result = first_circle.Intersect(second_circle);
		REQUIRE(result.size() == 0);
	}

}