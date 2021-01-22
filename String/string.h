#include <iostream>
#include <cstring>

class iterator {
 public:
  iterator() = default;

  iterator(int val, bool isnormal): value_(val), normal_(isnormal) {}

  ~iterator() = default;

  iterator& operator++ () {
    if (normal_) {
      ++value_;
    } else {
      --value_;
    }
    return *this;
  }

  bool operator==(const int& val) const {
    return value_ == val;
  }

  int value_number() const {
    return value_;
  }

 private:
  int value_{0};
  bool normal_{true};
};

class String {
 public:
  String() = default;

  String(const char* char_ptr) : size_(strlen(char_ptr)), string_(new char[2 * strlen(char_ptr)]),
                                          capacity_(2 * size_) {
    memcpy(string_, char_ptr, size_);
  }

  String(size_t size_value, char char_value = '\0') : size_(size_value), string_(new char[2 * size_value]),
                                                               capacity_(2 * size_) {
    memset(string_, char_value, size_value);
  }

  String(std::initializer_list<char> list) : size_(list.size()), string_(new char[2 * list.size()]),
                                             capacity_(2 * size_) {
    std::copy(list.begin(), list.end(), string_);
  }

  String(const String& input_string) : String(input_string.size_, '\0') {
    memcpy(string_, input_string.string_, size_);
  }

  void swap(String input_string) {
    std::swap(size_, input_string.size_);
    std::swap(string_, input_string.string_);
    std::swap(capacity_, input_string.capacity_);
  }

  String& operator=(String input_string) {
    swap(input_string);
    return* this;
  }

  char& operator[](size_t index) {
    return string_[index];
  }

  const char& operator[](size_t index) const {
    return string_[index];
  }

  size_t length() const {
    return size_;
  }

  void check_size_bigger() {
    if (size_ > capacity_) {
      capacity_ *= 2;
      char* copy_char = new char[capacity_];
      memcpy(copy_char, string_, size_);
      delete[] string_;
      string_ = copy_char;
    }
  }

  void push_back(const char char_input) {
    ++size_;
    check_size_bigger();
    string_[size_ - 1] = char_input;
  }

  void check_size_smaller() {
    if (capacity_ > size_ * 4) {
      capacity_ /= 2;
      char* copy_char = new char[capacity_];
      memcpy(copy_char, string_, size_);
      delete[] string_;
      string_ = copy_char;
    }
  }

  void pop_back() {
    --size_;
    check_size_smaller();
  }

  char& front() {
    return string_[0];
  }

  const char& front() const {
    return string_[0];
  }

  char& back() {
    if (size_ == 0) {
      return string_[0];
    } else {
      return string_[size_ - 1];
    }
  }

  const char& back() const {
    if (size_ == 0) {
      return string_[0];
    } else {
      return string_[size_ - 1];
    }
  }

  String& operator+=(const char input_char) {
    push_back(input_char);
    return *this;
  }

  String& operator+=(const String& input_string) {
    const size_t size_old = size_;
    size_ = size_ + input_string.length();
    check_size_bigger();
    memcpy(string_ + size_old, input_string.string_, input_string.length());
    return *this;
  }

  size_t find(const String& substring) const {
    return find_universal(substring, true);
  }

  size_t rfind(const String& substring) const {
    return find_universal(substring, false);
  }


  String substr(const size_t index, const size_t count) const {
    String copy_string;
    copy_string.size_ = count;
    copy_string.capacity_ = count * 2;
    copy_string.string_ = new char[capacity_];
    memcpy(copy_string.string_, string_ + index, count);
    return copy_string;
  }

  bool empty() const {
    return size_ == 0;
  }

  void clear() {
    size_ = 0;
    capacity_ = 1;
    delete [] string_;
    string_ = new char[1];
  }

  bool operator==(const String& input_string) const {
    bool ans = size_ == input_string.length();
    for (size_t i = 0; i < size_ && ans; i++) {
      if (string_[i] != input_string[i]) {
        ans = false;
      }
    }
    return ans;
  }

  ~String() {
    delete [] string_;
  }

 private:
  size_t size_{0};
  char* string_{new char[1]};
  size_t capacity_{1};

  size_t find_universal(const String& substring, bool is_normal) const { //is_normal == true for find, false otherwise
    if (size_ < substring.length()) {
      return size_;
    }
    bool found = false;
    int start_value;
    int end_value;
    if (is_normal) {
      start_value = 0;
      end_value = static_cast<int>(size_ + 1 - substring.length());
    } else {
      start_value = static_cast<int>(size_ - substring.length());
      end_value = -1;
    }
    iterator ans(start_value, is_normal);
    while (!(ans == end_value)) {
      bool found_local = true;
      size_t i = ans.value_number();
      while (i < static_cast<size_t>(ans.value_number() + substring.length())) {
        if (string_[i] != substring[static_cast<size_t>(i - ans.value_number())]) {
          found_local = false;
          break;
        }
        ++i;
      }
      found = found_local;
      if (found) {
        break;
      }
      ++ans;
    }
    if (found) {
      return static_cast<size_t>(ans.value_number());
    } else {
      return size_;
    }
  }
};

String operator+(const String& input_string, const char input_char) {
  String copy_string = input_string;
  copy_string += input_char;
  return copy_string;
}

String operator+(const char input_char, const String& input_string) {
  String copy_string(1, input_char);
  copy_string += input_string;
  return copy_string;
}

String operator+(const String& input_string_1, const String& input_string_2) {
  String copy_string = input_string_1;
  copy_string += input_string_2;
  return copy_string;
}

std::ostream& operator<<(std::ostream& out, const String& str) {
  for (size_t i = 0; i < str.length(); ++i) {
    out << str[i];
  }
  return out;
}

std::istream& operator>>(std::istream& in, String& str) {
  char char_input;
  str.clear();
  in >> std::noskipws;
  while (in >> char_input && char_input != '\n') {
    str.push_back(char_input);
  }
  return in;
}

