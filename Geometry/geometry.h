#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>
#include <cassert>

const double PI = 3.14159265358979323846264338327950288419716939937510;

class Line;

bool double_compare(double first, double second) {
  return (std::abs(first - second) < 0.000001);
}

struct Point{
 public:
  double x{0};
  double y{0};

  ~Point() = default;

  Point() = default;

  Point(const double& first, const double& second): x(first), y(second) {}

  Point(const Line& first, const Line& second);

  Point& operator=(Point val) {
    std::swap(x, val.x);
    std::swap(y, val.y);
    return *this;
  }

};

bool operator==(const Point& first, const Point& second) {
  return (double_compare(first.x, second.x) && double_compare(first.y,second.y));
}

bool operator!=(const Point& first, const Point& second) {
  return !(first == second);
}

double dist(const Point& a, const Point& b) {
  return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

class Vector {
 public:
  ~Vector() = default;

  Vector() = default;

  Vector(const Point& val) : vector_point(val) {}

  Vector& operator=(Vector vec) {
    std::swap(vector_point.x, vec.vector_point.x);
    std::swap(vector_point.y, vec.vector_point.y);
    return *this;
  }

  void rotate(double angle) {
    angle *= PI / 180.0;
    Point copy(vector_point.x * cos(angle) - vector_point.y * sin(angle), vector_point.x * sin(angle) + vector_point.y * cos(angle));
    vector_point = copy;
  }

  void rotatep(Point x, double angle) {
    Vector t = *this - x;
    t.rotate(angle);
    *this = x + t;
  }

  Vector& operator*=(const double& val) {
    vector_point.x *= val;
    vector_point.y *= val;
    return *this;
  }

  friend Vector operator+(Vector first, Vector second);

  friend Vector operator-(Vector first, Vector second);

  friend bool operator==(Vector first, Vector second);

  explicit operator Point() {
    return vector_point;
  }

  void reflex(Point center) {
    Vector t = *this - center;
    *this = center - t;
  }

  void reflex(Line& axis);

  double dist(Line& axis) const;

  void scale(Point center, double coefficient) {
    Vector t = *this - Vector(center);
    t *= coefficient;
    *this = Vector(center) + t;
  }

 private:
  Point vector_point{Point(0, 0)};
};

Vector operator+ (Vector first, Vector second) {
  return Vector(Point(first.vector_point.x + second.vector_point.x, first.vector_point.y + second.vector_point.y));
}

Vector operator- (Vector first, Vector second) {
  return Vector(Point(first.vector_point.x - second.vector_point.x, first.vector_point.y - second.vector_point.y));
}

bool operator==(Vector first, Vector second) {
  return (Point(first) == Point(second));
}

bool operator!=(Vector first, Vector second) {
  return (Point(first) != Point(second));
}


class Line{
 public:
  friend Vector;
  friend Point;

  ~Line() = default;

  Line(): a(1), b(0), c(0) {}

  Line(const double& a_val, const double& b_val, const double& c_val): a(a_val), b(b_val), c(c_val) {}

  Line(const Point& first, const Point& second) {
    a = first.y - second.y;
    b = second.x - first.x;
    c = first.x * second.y - second.x * first.y;
  }

  Line(const double& angle, const double& shift): a(angle), b(-1), c(shift) {}

  Line(const Point& val, const double& angle): a(angle), b(-1), c(val.y - val.x * angle) {}

  Line(const Line& l, const Point& p) {
    a = -l.b;
    b = l.a;
    c = p.x * l.b - p.y * l.a;
  }

  int side(const Point& value) const { // -1 (beneath the line), 0 (on the line) and 1 (above the line)
    double t = value.x * a + value.y * b + c;
    if (double_compare(t, 0)) {
      return 0;
    }
    if (t > 0) {
      return 1;
    } else {
      return -1;
    }
  }
  //line is defined by general equation: ax + by + c = 0
 private:
  double a;
  double b;
  double c;
  friend bool operator==(const Line& first, const Line& second);
};

bool operator==(const Line& first, const Line& second) {
  bool condition_one = double_compare(first.a * second.b, first.b * second.a);
  bool condition_two = double_compare(first.b * second.c, first.c * second.b);
  if (condition_one && condition_two) {
    return true;
  } else {
    return false;
  }
}

bool operator!=(const Line& first, const Line& second) {
  return !(first == second);
}

void Vector::reflex(Line& axis) {
  double temp = -2 * (axis.a * vector_point.x + axis.b * vector_point.y + axis.c) / (axis.a * axis.a + axis.b * axis.b);
  vector_point.x += temp * axis.a;
  vector_point.y += temp * axis.b;
}

double Vector::dist(Line& axis) const { //returns distance between vector_point and axis
  double ret = (axis.a * vector_point.x + axis.b * vector_point.y + axis.c) / sqrt((axis.a * axis.a + axis.b * axis.b));
  if (ret < 0) {
    ret *= -1;
  }
  return ret;
}

Point::Point(const Line& first, const Line& second) { //Intersection point of first and second lines
  double determinant = first.a * second.b - first.b * second.a;
  x = (first.b * second.c - first.c * second.b) / determinant;
  y = (first.c * second.a - first.a * second.c) / determinant;
  if (double_compare(x, 0)) {
    x = 0;
  }
  if (double_compare(y, 0)) {
    y = 0;
  }
}

class Shape {
 public:
  virtual ~Shape() = 0;

  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool isCongruentTo(const Shape& another) const = 0;

  virtual bool isSimilarTo(const Shape& another) const = 0;

  virtual bool containsPoint(Point point) const = 0;

  virtual void rotate(Point center, double angle) = 0;

  virtual void reflex(Point center) = 0;

  virtual void reflex(Line axis) = 0;

  virtual void scale(Point center, double coefficient) = 0;

  virtual bool operator== (const Shape& x) const = 0;

};

Shape::~Shape() {}

class Polygon: public Shape {
 public:
  ~Polygon() override = default;

  Polygon() = default;

  explicit Polygon(const std::vector<Point>& value): vertices_(value) {}

  Polygon(std::initializer_list<Point> list) {
    vertices_.resize(list.size());
    size_t i = 0;
    for (auto x : list) {
      vertices_[i] = x;
      ++i;
    }
  }

  size_t verticesCount() const {
    return vertices_.size();
  }

  const std::vector<Point>& getVertices() const {
    return vertices_;
  }

  bool isConvex() const {
    bool ret_val = true;
    for (size_t i = 0; i < verticesCount() && ret_val; ++i) {
      Point a = vertices_[i];
      Point b = vertices_[(i + 1) % verticesCount()];
      Point c = vertices_[(i + 2) % verticesCount()];
      Point d = vertices_[(i + 3) % verticesCount()];
      Line temp(b, c);
      ret_val = temp.side(a) == temp.side(d);
    }
    return ret_val;
  }

  virtual double perimeter() const override {
    double ans = 0;
    for (size_t i = 0; i < verticesCount() - 1; ++i) {
      ans += dist(vertices_[i], vertices_[i + 1]);
    }
    ans += dist(vertices_[0], vertices_[verticesCount() - 1]);
    return ans;
  }

  virtual double area() const override {
    double ans = 0;
    for (size_t i = 0; i < verticesCount(); ++i) {
      ans += vertices_[i].x * vertices_[(i + 1) % verticesCount()].y;
    }
    for (size_t i = 0; i < verticesCount(); ++i) {
      ans -= vertices_[i].y * vertices_[(i + 1) % verticesCount()].x;
    }
    ans /= 2;
    if (ans < 0) {
      ans *= -1;
    }
    return ans;
  }

  bool operator==(const Shape& another) const override {
    const auto* copy = dynamic_cast<const Polygon*>(&another);
    if (copy == nullptr) {
      return false;
    }
    if (copy->verticesCount() != verticesCount()) {
      return false;
    }
    size_t n = verticesCount();
    for (size_t i = 0; i < n; ++i) {
      bool f2 = true;
      for (size_t j = 0; j < n && f2; ++j) {
        if (vertices_[j] != copy->vertices_[(j + i) % n]) {
          f2 = false;
        }
      }
      if (f2) {
        return true;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      bool f2 = true;
      for (size_t j = 0; j < n && f2; ++j) {
        if (vertices_[j] != copy->vertices_[((n - j) + i) % n]) {
          f2 = false;
        }
      }
      if (f2) {
        return true;
      }
    }
    //здеаь ловит WA
    return false;
  }

  bool operator!=(const Shape& another) const {
    return !(*this == another);
  }

  bool isCongruentTo(const Shape& another) const override {
    const auto* copy = dynamic_cast<const Polygon*>(&another);

    if (copy == nullptr) {
      return false;
    }
    if (copy->verticesCount() != verticesCount()) {
      return false;
    }
    if (!double_compare(copy->perimeter(), perimeter())) {
      return false;
    }
    if (isSimilarTo(*copy)) {
      return true;
    }
    return false;
  }

  bool isSimilarTo(const Shape& another) const override {
    const auto* copy = dynamic_cast<const Polygon*>(&another);
    if (copy == nullptr) {
      return false;
    }
    if (copy->verticesCount() != verticesCount()) {
      return false;
    }
    size_t n = verticesCount();
    double similarity_rate;
    for (size_t i = 0; i < n; ++i) {
      bool local_ans = true;
      similarity_rate = 0;
      for (size_t j = 0; j < n; ++j) {
        size_t cur = (i + j);
        double a = dist(vertices_[cur % n], vertices_[(cur + 1) % n]);
        double b = dist(vertices_[(cur + 1) % n], vertices_[(cur + 2) % n]);
        double a1 = dist(copy->vertices_[j], copy->vertices_[(j + 1) % n]);
        double b1 = dist(copy->vertices_[(j + 1) % n], copy->vertices_[(j + 2) % n]);
        Line l(vertices_[cur], vertices_[(cur + 1) % n]);
        double h = Vector(vertices_[(cur + 2) % n]).dist(l);
        Line l1(copy->vertices_[j], copy->vertices_[(j + 1) % n]);
        double h1 = Vector(copy->vertices_[(j + 2) % n]).dist(l1);
        if (similarity_rate == 0) {
          similarity_rate = a / a1;
        }
        if (!double_compare(a * b1, b * a1) || !double_compare(h * a1, h1 * a) || !double_compare(a, a1 * similarity_rate)) {
          local_ans = false;
          break;
        }
      }
      if (local_ans) {
        return true;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      bool local_ans = true;
      similarity_rate = 0;
      for (size_t j = 0; j < n; ++j) {
        size_t cur = (n + i - j);
        double b = dist(vertices_[cur % n], vertices_[(cur + 1) % n]);
        double a = dist(vertices_[(cur + 1) % n], vertices_[(cur + 2) % n]);
        double a1 = dist(copy->vertices_[j], copy->vertices_[(j + 1) % n]);
        double b1 = dist(copy->vertices_[(j + 1) % n], copy->vertices_[(j + 2) % n]);
        Line l(vertices_[(cur + 2) % n], vertices_[(cur + 1) % n]);
        double h = Vector(vertices_[cur % n]).dist(l);
        Line l1(copy->vertices_[j], copy->vertices_[(j + 1) % n]);
        double h1 = Vector(copy->vertices_[(j + 2) % n]).dist(l1);
        if (similarity_rate == 0) {
          similarity_rate = a / a1;
        }
        if (!double_compare(a * b1, b * a1) || !double_compare(h * a1, h1 * a) ||
            !double_compare(a, a1 * similarity_rate)) {
          local_ans = false;
          break;
        }
      }
      if (local_ans) {
        return true;
      }
    }
    return false;
  }

  bool containsPoint(Point point) const override {
    size_t n = verticesCount();
    size_t chet = 0;
    for (size_t i = 0; i < n; ++i) {
      Point a = vertices_[(i + n - 1) % n];
      Point b = vertices_[i];
      bool f1 = point.x > Point(Line(a, b), Line(0, -1, point.y)).x;
      bool f2 = (point.y > std::min(a.y, b.y) && point.y < std::max(a.y, b.y));
      if (f1 && f2) {
        ++chet;
      }
    }
    return (chet % 2 == 1);
  }

  void rotate(Point center, double angle) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      Vector t(vertices_[i]);
      t.rotatep(center, angle);
      vertices_[i] = Point(t);
    }
  }

  void reflex(Point center) override {
    for (size_t i = 0; i< verticesCount(); ++i) {
      Vector t(vertices_[i]);
      t.reflex(center);
      vertices_[i] = Point(t);
    }
  }

  void reflex(Line axis) override {
    for (size_t i = 0; i< verticesCount(); ++i) {
      Vector t(vertices_[i]);
      t.reflex(axis);
      vertices_[i] = Point(t);
    }
  }

  void scale(Point center, double coefficient) override {
    for (size_t i = 0; i< verticesCount(); ++i) {
      Vector t(vertices_[i]);
      t.scale(center, coefficient);
      vertices_[i] = Point(t);
    }
  }

 protected:
  std::vector<Point> vertices_;
};

class Ellipse : public Shape {
 public:
  ~Ellipse() override = default;

  Ellipse() = default;

  Ellipse(const Point& focus_one, const Point& focus_two, const double& dist): a(dist / 2), f1(focus_one), f2(focus_two) {}

  std::pair<Point, Point> focuses() const{
    return {f1, f2};
  }

  std::pair<Line, Line> directrises() const{
    double c = dist(f1, f2) / 2;
    double t = a * a / c;
    Point cent = center();
    if (f1.y == f2.y) {
      return {Line(1, 0, -cent.x - t), Line(1, 0, -cent.x + t)};
    }
    if (f1.x == f2.x) {
      return {Line(0, 1, -cent.y - t), Line(0, 1, -cent.y + t)};
    }
    double k = (f1.y - f2.y) / (f1.x - f2.x);
    double sink = 1.0 / (1.0 + k * k);
    double cosk = 1 - sink;
    sink = sqrt(sink);
    cosk = sqrt(cosk);
    if (k < 0) {
      cosk *= -1;
    }
    Point d1(cent.x + cosk * t, cent.y + sink * t);
    Point d2(cent.x - cosk * t, cent.y - sink * t);
    return {Line(d1, -1.0/k), Line(d2, -1.0/k)};
  }

  double eccentricity() const {
    double c = dist(f1, f2) / 2;
    return c / a;
  }

  Point center() const {
    return Point((f1.x + f2.x) / 2, (f1.y + f2.y) / 2);
  }

  virtual double perimeter() const override {
    return 4 * a * std::comp_ellint_2(eccentricity());
  }

  virtual double area() const override {
    double c = dist(f1, f2) / 2;
    double b = sqrt(a * a - c * c);
    return PI * a * b;
  }

  bool operator==(const Shape& another) const override {
    const auto* copy = dynamic_cast<const Ellipse*>(&another);
    if (copy == nullptr) {
      return false;
    }
    if (copy->a == a && ((copy->f1 == f1 && copy->f2 == f2) || (copy->f1 == f2 && copy->f2 == f1))) {
      return true;
    } else {
      return false;
    }
  }

  bool isCongruentTo(const Shape& another) const override {
    const auto* copy = dynamic_cast<const Ellipse*>(&another);
    if (copy == nullptr) {
      return false;
    }
    double c = dist(f1, f2) / 2;
    double c1 = dist(copy->f1, copy->f2) / 2;
    double b = sqrt(a * a - c * c);
    double b1 = sqrt(copy->a * copy-> a - c1 * c1);
    if (double_compare(copy->a, a) && double_compare(b, b1)) {
      return true;
    }
    if (double_compare(copy->a, b) && double_compare(a, b1)) {
      return true;
    }
    return false;
  }

  bool isSimilarTo(const Shape& another) const override {
    const auto* copy = dynamic_cast<const Ellipse*>(&another);
    if (copy == nullptr) {
      return false;
    }
    double c = dist(f1, f2) / 2;
    double c1 = dist(copy->f1, copy->f2) / 2;
    double b = sqrt(a * a - c * c);
    double b1 = sqrt(copy->a * copy-> a - c1 * c1);
    if (double_compare(b * copy->a, b1 * a)) {
      return true;
    }
    if (double_compare(a * copy->a, b * b1)) {
      return true;
    }
    return false;
  }

  bool containsPoint(Point point) const override {
    return (dist(point, f1) + dist(point, f2) < 2 * a);
  }

  void rotate(Point center, double angle) override {
    Vector t(f1);
    t.rotatep(center, angle);
    f1 = Point(t);
    t = Vector(f2);
    t.rotatep(center, angle);
    f2 = Point(t);
  }

  void reflex(Point center) override {
    Vector t(f1);
    t.reflex(center);
    f1 = Point(t);
    t = Vector(f2);
    t.reflex(center);
    f2 = Point(t);
  }

  void reflex(Line axis) override {
    Vector t(f1);
    t.reflex(axis);
    f1 = Point(t);
    t = Vector(f2);
    t.reflex(axis);
    f2 = Point(t);
  }

  void scale(Point center, double coefficient) override {
    a *= coefficient;
    if (a < 0) {
      a *= -1;
    }
    Vector t(f1);
    t.scale(center, coefficient);
    f1 = Point(t);
    t = Point(f2);
    t.scale(center, coefficient);
    f2 = Point(t);
  }

 protected:
  double a{1}; //semi-major axis

 private:
  Point f1{Point(0, 0)}; // first_focus
  Point f2{Point(0, 0)}; // second_focus
};

class Circle: public Ellipse {
 public:
  ~Circle() override = default;

  Circle(): Ellipse(Point(0,0), Point(0,0), 0.5) {}

  Circle(const Point& val, const double& rad): Ellipse(val, val, 2*rad) {}

  double radius() const {
    return a;
  }
};

class Rectangle: public Polygon {
 public:
  ~Rectangle() override = default;

  using Polygon::Polygon;

  Rectangle(const Point& first, const Point& second, double ratio) {
    if (ratio < 1) {
      ratio = 1 / ratio;
    }
    Vector diagonal = Vector(second) - Vector(first);
    diagonal.rotate(atan(ratio) * 180 / PI);
    diagonal *= 1.0 / sqrt(1.0 + ratio * ratio);
    vertices_.resize(4);
    vertices_[0] = first;
    vertices_[1] = Point(Vector(first) + diagonal);
    vertices_[2] = second;
    diagonal *= -1;
    vertices_[3] = Point(Vector(second) + diagonal);
  }

  Point center() const {
    Point ret_val;
    ret_val.x = (vertices_[0].x + vertices_[1].x + vertices_[2].x + vertices_[3].x) / 4;
    ret_val.y = (vertices_[0].y + vertices_[1].y + vertices_[2].y + vertices_[3].y) / 4;
    return ret_val;
  }

  std::pair<Line, Line> diagonals() const {
    return {Line(vertices_[0], vertices_[2]), Line(vertices_[1], vertices_[3])};
  }
};

class Square: public Rectangle {
 public:
  ~Square() override = default;

  using Rectangle::Rectangle;

  Square(const Point& first, const Point& second): Rectangle(first, second, 1.0) {}

  Circle circumscribedCircle()  const {
    double r = dist(vertices_[0], vertices_[2]) / 2;
    return Circle(center(), r);
  }

  Circle inscribedCircle() const {
    double r = dist(vertices_[0], vertices_[1]) / 2;
    return Circle(center(), r);
  }
};

class Triangle: public Polygon {
 public:
  ~Triangle() override = default;

  Triangle(Point a, Point b, Point c) {
    vertices_.push_back(a);
    vertices_.push_back(b);
    vertices_.push_back(c);
  }

  using Polygon::Polygon;

  Circle circumscribedCircle() const {
    Point a = vertices_[0];
    Point b = vertices_[1];
    Point c = vertices_[2];
    Line ab(a, b);
    Line bc(b, c);
    Point cent;
    Line s1(ab, Point((Vector(a) + Vector(b)) *= 0.5));
    Line s2(bc, Point((Vector(b) + Vector(c)) *= 0.5));
    cent = Point(s1, s2);
    double radi = dist(cent, a);
    return Circle(cent, radi);
  }

  Circle inscribedCircle() const {
    Point a = vertices_[0];
    Point b = vertices_[1];
    Point c = vertices_[2];
    double ab = dist(a, b);
    double bc = dist(b, c);
    double ac = dist(a, c);
    Point cent;
    Line bi_1(a, Point(((Vector(b) *= ac) + (Vector(c) *= ab)) *= (1 / (ac + ab))));
    Line bi_2(b, Point(((Vector(a) *= bc) + (Vector(c) *= ab)) *= (1 / (bc + ab))));
    cent = Point(bi_1, bi_2);
    Line temp = Line(a, b);
    double radi = Vector(cent).dist(temp);
    return Circle(cent, radi);
  }

  Point centroid() const {
    Point a = vertices_[0];
    Point b = vertices_[1];
    Point c = vertices_[2];
    return Point((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3);
  }

  Point orthocenter() const {
    Point a = vertices_[0];
    Point b = vertices_[1];
    Point c = vertices_[2];
    Point ans;
    ans.x = (b.x * (a.x - c.x) + b.y * (a.y - c.y)) * (c.y - b.y) - (c.y - a.y) * (a.x * (b.x - c.x) + a.y * (b.y - c.y));
    ans.x /= (c.x - b.x) * (c.y - a.y) - (c.y - b.y) * (c.x - a.x);
    ans.y = (b.x * (a.x - c.x) + b.y * (a.y - c.y)) * (c.x - b.x) - (c.x - a.x) * (a.x * (b.x - c.x) + a.y * (b.y - c.y));
    ans.y /= (c.y - b.y) * (c.x - a.x) - (c.x - b.x) * (c.y - a.y);
    return ans;
  }

  Line EulerLine() const {
    return Line(orthocenter(), circumscribedCircle().center());
  }

  Circle ninePointsCircle() const {
    Point a = vertices_[0];
    Point b = vertices_[1];
    Point c = vertices_[2];
    Point a1((a.x + b.x) / 2, (a.y + b.y) / 2);
    Point b1((b.x + c.x) / 2, (b.y + c.y) / 2);
    Point c1((c.x + a.x) / 2, (c.y + a.y) / 2);
    return Triangle({a1, b1, c1}).circumscribedCircle();
  }
};


bool operator!= (const Shape& a, const Shape& b) {
  return !(a == b);
}
