#include <iostream>
#include <vector>
#include <string>

class BigInteger {
 public:

  ~BigInteger() = default;

  BigInteger() : BigInteger(0) {}

  BigInteger(long value) {
    if (value >= 0) {
      positive_ = true;
    } else {
      value *= -1;
      positive_ = false;
    }
    while (value > 0) {
      numbers_.push_back(value % base_);
      value /= base_;
    }
    if (numbers_.empty()) {
      numbers_.push_back(0);
    }
  }

  BigInteger(const BigInteger& value) : numbers_(value.numbers_), positive_(value.positive_) {}

  BigInteger(const std::string& value) { //+
    long starting_digit = 0;
    if (value[0] == '-') {
      starting_digit = 1;
      positive_ = false;
    } else {
      starting_digit = 0;
      positive_ = true;
    }
    long long multiplier = 1;
    long digit_number = 0;
    long long dig = 0;
    numbers_.resize((value.size() - starting_digit + base_number_ - 1) / base_number_);
    for (long i = static_cast<long>(value.size()) - 1; i >= starting_digit; --i) {
      dig += (value[i] - '0') * multiplier;
      multiplier *= 10;
      if (multiplier == base_ || i == starting_digit) {
        multiplier = 1;
        numbers_[digit_number] = dig;
        dig = 0;
        ++digit_number;
      }
    }
    normalize();
  }

  BigInteger& operator=(const long value) {
    *this = std::to_string(value);
    return *this;
  }

  BigInteger& operator=(const std::string& value) {
    *this = BigInteger(value);
    return *this;
  }

  void swap(BigInteger val) {
    std::swap(numbers_, val.numbers_);
    std::swap(positive_, val.positive_);
  }

  BigInteger& operator=(const BigInteger& val) {
    swap(val);
    return *this;
  }

  void normalize() {
    long i = static_cast<long>(numbers_.size()) - 1;
    while (i > 0 && numbers_[i] == 0) {
      --i;
      numbers_.pop_back();
    }
    if (numbers_.size() == 1 && numbers_[0] == 0) {
      positive_ = true;
    }
  }

  BigInteger operator-() const;

  BigInteger& operator++() {
    *this += 1;
    return *this;
  }

  BigInteger& operator--() {
    *this -= 1;
    return *this;
  }

  BigInteger operator++(int) {
    BigInteger copy = *this;
    ++*this;
    return copy;
  }

  BigInteger operator--(int) {
    BigInteger copy = *this;
    --*this;
    return copy;
  }

  BigInteger& operator+=(long value) {
    if (value < 0) {
      *this -= -value;
      return *this;
    }
    long carry = value;
    long temp;
    for (long i = 0; i < static_cast<long>(numbers_.size()) && value > 0; ++i) {
      temp = (carry + numbers_[i]) % base_;
      carry = (carry + numbers_[i]) / base_;
      numbers_[i] = temp;
    }
    if (carry > 0) {
      numbers_.push_back(carry);
    }
    normalize();
    return *this;
  }

  BigInteger& operator+=(const BigInteger& value);

  BigInteger& operator-=(const long value) {
    BigInteger temp(value);
    *this -= temp;
    return *this;
  }

  BigInteger& operator-=(const BigInteger& value);

  BigInteger& operator*=(const int value) {
    if (value == 0) {
      *this = 0;
      return *this;
    }
    if ((value >= 0) == positive_) {
      positive_ = true;
    } else {
      positive_ = false;
    }
    long carry = 0;
    long val = std::abs(value);
    for (long long& number : numbers_) {
      long copy = val * static_cast<long>(number);
      copy += carry;
      carry = copy / static_cast<long>(base_);
      number = copy % static_cast<long>(base_);
    }
    if (carry > 0) {
      numbers_.push_back(carry);
    }
    normalize();
    return *this;
  }

  void add_nulls(const long null_amount) {
    if (null_amount > 0) {
      long t = numbers_.size();
      numbers_.resize(t + null_amount);
      for (long i = static_cast<long>(numbers_.size()) - 1; i > -1; --i) {
        numbers_[i] = numbers_[i - null_amount];
      }
      for (long i = 0; i < null_amount; ++i) {
        numbers_[i] = 0;
      }
    }
  }

  void del_last() {
    for (size_t i = 0; i < numbers_.size() - 1; ++i) {
      numbers_[i] = numbers_[i + 1];
    }
    numbers_.pop_back();
  }

  BigInteger& operator*=(const BigInteger& value) {
    bool positiv;
    if (positive_ != value.positive_) {
      positiv = false;
    } else {
      positiv = true;
    }
    BigInteger copy;
    BigInteger sum = 0;
    for (size_t i = 0; i < value.numbers_.size(); ++i) {
      copy = *this;
      copy *= value.numbers_[i];
      copy.add_nulls(i);
      copy.normalize();
      sum += copy;
    }
    *this = sum;
    this->positive_ = positiv;
    normalize();
    return *this;
  }

  BigInteger abs() const {
    BigInteger copy = *this;
    copy.positive_ = true;
    return copy;
  }

  BigInteger& operator/=(const BigInteger& value);

  BigInteger& operator%=(const int value) {
    BigInteger copy = value;
    *this %= copy;
    return *this;
  }

  BigInteger& operator%=(const BigInteger& value);

  void Clear() {
    numbers_.clear();
    positive_ = true;
  }

  std::string toString() const;

  explicit operator bool() const;

  explicit operator double() const {
    double res = 0;
    double multiplier = 1;
    for (long digit : numbers_) {
      res += digit * multiplier;
      multiplier *= base_;
    }
    if (!positive_) {
      res *= -1;
    }
    return res;
  }

 private:
  std::vector<long long> numbers_;
  bool positive_{true};
  const long long base_ = 10;
  const int base_number_ = 1;

  friend BigInteger operator+(const BigInteger& number_first, const BigInteger& number_second);

  friend BigInteger operator*(const BigInteger& number_first, const BigInteger& number_second);

  friend std::ostream& operator<<(std::ostream& out, const BigInteger& value);

  friend std::ostream& operator<<(std::ostream& out, const BigInteger& value);

  friend bool operator<(const BigInteger& first_number, const BigInteger& second_number);

  friend class Rational;

  std::pair<BigInteger, BigInteger> div_and_mod(const BigInteger& value);
};


bool operator<(const BigInteger& first_number, const BigInteger& second_number) { //+
  if (first_number.positive_ && !second_number.positive_) {
    return false;
  }
  if (!first_number.positive_ && second_number.positive_) {
    return true;
  }
  bool same_sign;
  if (first_number.positive_ && second_number.positive_) {
    same_sign = true;
  } else {
    same_sign = false;
  }
  if (first_number.numbers_.size() < second_number.numbers_.size()) {
    return same_sign;
  }
  if (first_number.numbers_.size() > second_number.numbers_.size()) {
    return !same_sign;
  }
  if (first_number.numbers_.size() == second_number.numbers_.size()) {
    for (long i = static_cast<long>(first_number.numbers_.size()) - 1; i > -1; --i) {
      if (first_number.numbers_[i] > second_number.numbers_[i]) {
        return !same_sign;
      }
      if (first_number.numbers_[i] < second_number.numbers_[i]) {
        return same_sign;
      }
    }
    return false;
  }
  return true;
}

bool operator==(const BigInteger& first_number, const BigInteger& second_number) {//+
  return !((first_number < second_number) || (second_number < first_number));
}

bool operator>(const BigInteger& first_number, const BigInteger& second_number) {//+
  return (second_number < first_number);
}

bool operator<=(const BigInteger& first_number, const BigInteger& second_number) {//+
  return !(first_number > second_number);
}

bool operator>=(const BigInteger& first_number, const BigInteger& second_number) {//+
  return !(first_number < second_number);
}

bool operator!=(const BigInteger& first_number, const BigInteger& second_number) {//+
  return !(first_number == second_number);
}

BigInteger BigInteger::operator-() const { //+
  BigInteger ret_val = *this;
  if (ret_val != 0) {
    ret_val.positive_ = !positive_;
  }
  return ret_val;
}

BigInteger operator+(const BigInteger& number_first, const BigInteger& number_second) {
  BigInteger return_value;
  return_value = number_first;
  return_value += number_second;
  return return_value;
}

BigInteger operator-(const BigInteger& number_first, const BigInteger& number_second) {
  BigInteger return_value = number_first;
  return_value -= number_second;
  return return_value;
}

BigInteger operator*(const BigInteger& number_first, const BigInteger& number_second) {
  BigInteger return_value;
  return_value = number_first;
  return_value *= number_second;
  return return_value;
}

BigInteger operator/(const BigInteger& number_first, const BigInteger& number_second) {
  BigInteger return_value = number_first;
  return_value /= number_second;
  return return_value;
}

BigInteger operator%(const BigInteger& number_first, const BigInteger& number_second) {
  BigInteger return_value = number_first;
  return_value %= number_second;
  return return_value;
}


BigInteger& BigInteger::operator+=(const BigInteger& value) {//++
  if (value.positive_ != positive_) {
    if (positive_) {
      *this -= -value;
    } else {
      *this = value + *this;
    }
    return *this;
  }
  long carry = 0;
  numbers_.resize(std::max(numbers_.size(), value.numbers_.size()));
  size_t i;
  for (i = 0; i < std::min(numbers_.size(), value.numbers_.size()); ++i) {
    carry += numbers_[i] + value.numbers_[i];
    numbers_[i] = carry % base_;
    carry /= base_;
  }
  while (carry > 0) {
    if (i >= numbers_.size()) {
      numbers_.push_back(0);
    }
    carry += numbers_[i];
    numbers_[i] = carry % base_;
    carry /= base_;
    ++i;
  }
  normalize();
  return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& value) {//+
  if (!value.positive_) {
    *this += -value;
    return *this;
  }
  if (positive_) {
    if (*this >= value) {
      int carry = 0;
      size_t i = 0;
      for (i = 0; i < std::min(numbers_.size(), value.numbers_.size()); ++i) {
        if (numbers_[i] < value.numbers_[i] + carry) {
          numbers_[i] += base_;
          numbers_[i] -= value.numbers_[i];
          numbers_[i] -= carry;
          carry = 1;
        } else {
          numbers_[i] -= value.numbers_[i];
          numbers_[i] -= carry;
          carry = 0;
        }
      }
      while (carry > 0) {
        if (i >= numbers_.size()) {
          numbers_.push_back(0);
        }
        if (numbers_[i] < carry) {
          numbers_[i] += base_;
          numbers_[i] -= carry;
          carry = 1;
        } else {
          numbers_[i] -= carry;
          carry = 0;
        }
        ++i;
      }
    } else {
      *this = -(value - *this);
    }
  } else {
    positive_ = true;
    *this += value;
    positive_ = false;
  }
  normalize();
  return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& value) {
  *this = div_and_mod(value).first;
  return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& value) {
  *this = div_and_mod(value).second;
  return *this;
}

std::pair<BigInteger, BigInteger> BigInteger::div_and_mod(const BigInteger& value) {
  bool mod_positive = positive_;
  if (positive_ == value.positive_) {
    positive_ = true;
  } else {
    positive_ = false;
  }
  BigInteger copy = value;
  copy.positive_ = true;
  BigInteger res;
  res.positive_ = positive_;
  positive_ = true;
  if (*this < copy) {
    positive_ = mod_positive;
    return {0, *this};
  }
  res.numbers_.resize(numbers_.size() - value.numbers_.size() + 1);
  copy.add_nulls(static_cast<long>(numbers_.size() - value.numbers_.size()));
  for (long i = static_cast<long>(numbers_.size() - value.numbers_.size()); i > -1; --i) {
    long left = 0;
    long right = base_;
    while (right > left + 1) {
      long m = (right + left) / 2;
      if (*this >= copy * m) {
        left = m;
      } else {
        right = m;
      }
    }
    positive_ = mod_positive;
    res.numbers_[i] = left;
    *this -= left * copy;
    copy.del_last();
  }
  normalize();
  res.normalize();
  return {res, *this};
}

int digits_count(long value) {
  int ret_val = 0;
  if (value == 0) {
    return 1;
  }
  while (value > 0) {
    ++ret_val;
    value /= 10;
  }
  return ret_val;
}

std::string BigInteger::toString() const {
  std::string ret_val;
  if (numbers_.empty()) {
    return ret_val;
  }
  if (!positive_) {
    ret_val += "-";
  }
  ret_val += std::to_string(numbers_[numbers_.size() - 1]);
  for (long i = static_cast<long>(numbers_.size()) - 2; i > -1; --i) {
    for (int j = 0; j < base_number_ - digits_count(numbers_[i]); ++j) {
      ret_val += "0";
    }
    ret_val += std::to_string(numbers_[i]);
  }
  return ret_val;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& value) {
  out << value.toString();
  return out;
}

std::istream& operator>>(std::istream& in, BigInteger& value) {
  value.Clear();
  std::string input;
  in >> input;
  value = input;
  return in;
}

BigInteger::operator bool() const {
  return (*this != 0);
}

class Rational {
 public:
  ~Rational() = default;

  Rational() {
    positive_ = true;
  }

  Rational(const BigInteger& first, const BigInteger& second) {
    positive_ = first.positive_ == second.positive_;
    numerator_ = first;
    denominator_ = second;
    numerator_.positive_ = true;
    denominator_.positive_ = true;
    normalize();
  }

  Rational(const long first, const long second) {
    positive_ = (first >= 0) == (second >= 0);
    numerator_ = first;
    denominator_ = second;
    numerator_.positive_ = true;
    denominator_.positive_ = true;
    normalize();
  }

  Rational(const long value) {
    positive_ = value >= 0;
    numerator_ = std::abs(value);
    denominator_ = 1;
  }

  Rational(const BigInteger& value) {
    positive_ = value.positive_;
    numerator_ = value.abs();
    denominator_ = 1;
  }

  Rational& operator+=(const Rational& value) {
    long sign1, sign2;
    if (positive_) {
      sign1 = 1;
    } else {
      sign1 = -1;
    }
    if (value.positive_) {
      sign2 = 1;
    } else {
      sign2 = -1;
    }
    numerator_ = numerator_ * value.denominator_ * sign1 + denominator_ * value.numerator_ * sign2;
    denominator_ *= value.denominator_;
    normalize_sign();
    normalize();
    return *this;
  }

  Rational& operator-=(const Rational& value) {
    long sign1, sign2;
    if (positive_) {
      sign1 = 1;
    } else {
      sign1 = -1;
    }
    if (value.positive_) {
      sign2 = 1;
    } else {
      sign2 = -1;
    }
    numerator_ = numerator_ * value.denominator_ * sign1 - denominator_ * value.numerator_ * sign2;
    denominator_ *= value.denominator_;
    normalize_sign();
    normalize();
    return *this;
  }

  Rational& operator*=(const Rational& value) {
    long sign1, sign2;
    if (positive_) {
      sign1 = 1;
    } else {
      sign1 = -1;
    }
    if (value.positive_) {
      sign2 = 1;
    } else {
      sign2 = -1;
    }
    numerator_ *= (value.numerator_ * sign1 * sign2);
    denominator_ *= value.denominator_;
    normalize_sign();
    normalize();
    return *this;
  }

  Rational& operator/=(const Rational& value) {
    long sign1, sign2;
    if (positive_) {
      sign1 = 1;
    } else {
      sign1 = -1;
    }
    if (value.positive_) {
      sign2 = 1;
    } else {
      sign2 = -1;
    }
    numerator_ *= (value.denominator_ * sign1);
    denominator_ *= (value.numerator_ * sign2);
    normalize_sign();
    normalize();
    return *this;
  }

  friend bool operator<(const Rational& first, const Rational& second);

  Rational operator-() const {
    Rational copy = *this;
    copy.positive_ = !positive_;
    return copy;
  }

  std::string toString() const {
    std::string out;
    if (numerator_.numbers_.empty()) {
      return out;
    }
    if (!positive_) {
      out += "-";
    }
    out += numerator_.toString();
    if (denominator_ > 1) {
      out += "/";
      out += denominator_.toString();
    }
    return out;
  }

  std::string asDecimal(size_t precision = 0) {
    std::string out;
    if (numerator_.numbers_.empty()) {
      return out;
    }
    if (!positive_) {
      out += "-";
    }
    BigInteger copy1 = numerator_ / denominator_;
    if (precision == 0) {
      return copy1.toString();
    }
    out += copy1.toString();
    if (precision > 0) {
      out += ".";
      BigInteger copy2 = numerator_;
      copy2.add_nulls(precision);
      copy2.normalize();
      copy2 /= denominator_;
      copy1.add_nulls(precision);
      copy1.normalize();
      copy2 -= copy1;
      copy2.normalize();
      std::string str1 = copy2.toString();
      for (size_t i = 0; i < (precision - str1.length()); ++i) {
        out += "0";
      }
      out += str1;
    }
    return out;
  }

  explicit operator double() const {
    double res;
    res = static_cast<double>(numerator_) / static_cast<double>(denominator_);
    if (!positive_) {
      res *= -1;
    }
    return res;
  }

 private:
  BigInteger numerator_;
  BigInteger denominator_;
  bool positive_{true};

  void div_2(BigInteger& arg) {
    for (size_t i = 0; i < arg.numbers_.size(); ++i) {
      if (arg.numbers_[i] & 1) {
        if (i > 0) {
          arg.numbers_[i - 1] += arg.base_ / 2;
        }
      }
      arg.numbers_[i] /= 2;
    }
    for (size_t i = 0; i < arg.numbers_.size(); ++i) {
      if (arg.numbers_[i] > arg.base_) {
        arg.numbers_[i + 1] += arg.numbers_[i] / arg.base_;
        arg.numbers_[i] %= arg.base_;
      }
    }
    arg.normalize();
  }

  BigInteger GCD(BigInteger first_number, BigInteger second_number) { //first_number and second_number are positive!!!
    //this is code of binary gcd
    BigInteger multiplier = 1;
    while (true) {
      if (first_number == 0) {
        return (second_number * multiplier);
      }
      if (second_number == 0) {
        return (first_number * multiplier);
      }
      if (first_number == second_number) {
        return (first_number * multiplier);
      }
      if (first_number == 1) {
        return multiplier;
      }
      if (second_number == 1) {
        return multiplier;
      }
      if (second_number.numbers_[0] % 2 == 0 && first_number.numbers_[0] % 2 == 0) {
        div_2(first_number);
        div_2(second_number);
        multiplier *= 2;
        continue;
      }
      if (first_number.numbers_[0] % 2 == 0) {
        div_2(first_number);
        continue;
      }
      if (second_number.numbers_[0] % 2 == 0) {
        div_2(second_number);
        continue;
      }
      if (first_number > second_number) {
        first_number -= second_number;
        div_2(first_number);
        continue;
      }
      if (second_number > first_number) {
        second_number -= first_number;
        div_2(second_number);
        continue;
      }
    }
  }

  void normalize() {
    numerator_.positive_ = true;
    denominator_.positive_ = true;
    BigInteger del = GCD(numerator_, denominator_);
    numerator_ /= del;
    denominator_ /= del;
  }

  void normalize_sign() {
    positive_ = (numerator_.positive_ == denominator_.positive_) || (numerator_ == 0);
  }
};

Rational operator+(const Rational& first, const Rational& second) {
  Rational copy = first;
  copy += second;
  return copy;
}

Rational operator-(const Rational& first, const Rational& second) {
  Rational copy = first;
  copy -= second;
  return copy;
}

Rational operator*(const Rational& first, const Rational& second) {
  Rational copy = first;
  copy *= second;
  return copy;
}

Rational operator/(const Rational& first, const Rational& second) {
  Rational copy = first;
  copy /= second;
  return copy;
}

bool operator<(const Rational& first, const Rational& second) {
  if (first.positive_ && !second.positive_) {
    return false;
  }
  if (!first.positive_ && second.positive_) {
    return true;
  }
  bool same_sign;
  if (first.positive_) {
    same_sign = true;
  } else {
    same_sign = false;
  }
  if ((first.numerator_ * second.denominator_) < (first.denominator_ * second.numerator_)) {
    return same_sign;
  } else {
    if ((first.numerator_ * second.denominator_) > (first.denominator_ * second.numerator_)) {
      return !same_sign;
    } else {
      return false;
    }
  }
}

bool operator>=(const Rational& first, const Rational& second) {
  return !(first < second);
}

bool operator==(const Rational& first, const Rational& second) {
  return !(first < second || second < first);
}

bool operator!=(const Rational& first, const Rational& second) {
  return !(first == second);
}

bool operator>(const Rational& first, const Rational& second) {
  return (second < first);
}

bool operator<=(const Rational& first, const Rational& second) {
  return !(first > second);
}

