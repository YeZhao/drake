#include <iostream>
#include <string>
#include <vector>
#include <map>
//#include <Eigen/Core>

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {
///////////////////////////
// Function Declarations //
///////////////////////////

// JSON specific parsers
std::map<std::string, std::string> parse_json_object(const std::string&);
std::vector<std::string> parse_json_list(const std::string&);
template<int size>
Eigen::Matrix<double, size, 1>  json_to_eigen_d(const std::string&);
// Eigen::VectorXd json_to_eigen_i<int>(const std::string&);
double parse_json_double(const std::string&);
// misc string helpers
bool contains(const std::string&, const char*);
bool contains(const std::string&, char);
std::size_t find_closing_delim(const std::string&, std::size_t, char='{', char='}');


/////////////////////////
// Function Deinitions //
/////////////////////////
std::vector<std::string> parse_json_list(const std::string &data){
  std::vector<std::string> list;
  
  auto begin = data.find('[');
  begin = data.find_first_not_of(" \t\r\n",begin+1);
  auto end = begin;
  auto stop = find_closing_delim(data, begin,'[',']');

  while(begin < stop && begin != std::string::npos){
    if (data.at(begin) == '\"'){
      end = find_closing_delim(data,begin,'\"','\"');
      list.push_back(data.substr(begin+1,end-begin-1));
    }else if (data.at(begin) == '{'){
      end = find_closing_delim(data,begin,'{','}');
      list.push_back(data.substr(begin,end-begin+1));
    }else if (data.at(begin) == '['){
      end = find_closing_delim(data,begin,'[',']');
      list.push_back(data.substr(begin,end-begin+1));
    }else if(begin==data.find("Infinity", begin)){// special numerics
      list.push_back("Infinity");
      end = begin+7;
    }else if(begin==data.find("-Infinity", begin)){
      list.push_back("-Infinity");
      end = begin+8;
    }else if(begin==data.find("NaN",begin)){
      list.push_back("NaN");
      end = begin+2;
    }else if (contains("-0123456789.",data[begin])){// regular numerics
      end = data.find_first_not_of("-0123456789.e",begin)-1;// allow negative numbers, decimals and scientific notation
      list.push_back(data.substr(begin,end-begin+1));
    }else if (begin==data.find("null", begin)){
      list.push_back("null");
      end = begin+3;
    }else if (begin==data.find("true", begin)){
      list.push_back("true");
      end = begin+3;
    }else if (begin==data.find("false", begin)){
      list.push_back("false");
      end = begin+4;
    }else{
      // TODO define custom exception types
      // bad JSON format: the type cannot be recognized
      throw(1);
    }
    begin = data.find_first_not_of(", \t\n\r", end+1);
  }

  return list;
}

std::map<std::string, std::string> parse_json_object(const std::string &data){
  // initialize variables
  std::map<std::string, std::string> map;
  std::string key, val;
  std::size_t begin_key = data.find('{');
  std::size_t end_key, begin_val, end_val;
  
  while (begin_key != std::string::npos){

    // keys are all strings
    begin_key = data.find('\"',begin_key);
    if (begin_key == std::string::npos)
      break;
    begin_key++;
    end_key = data.find('\"',begin_key)-1;
    key = data.substr(begin_key,end_key-begin_key+1);
  
    // determine the value type
    begin_val = data.find_first_not_of(": \t\r\n",end_key+2);
    if (data.at(begin_val) == '\"'){
      end_val = find_closing_delim(data,begin_val,'\"','\"');
      begin_key = end_val+1;
      begin_val++; end_val--;
    }else if (data.at(begin_val) == '{'){
      end_val = find_closing_delim(data,begin_val,'{','}');
      begin_key = end_val+1;
    }else if (data.at(begin_val) == '['){
      end_val = find_closing_delim(data,begin_val,'[',']');
      begin_key = end_val+1;
    }else if(begin_val==data.find("Infinity", begin_val)){// special numerics
      end_val = begin_val+7;
      begin_key = end_val+1;
    }else if(begin_val==data.find("-Infinity", begin_val)){
      end_val = begin_val+8;
      begin_key = end_val+1;
    }else if(begin_val==data.find("NaN",begin_val)){
      end_val = begin_val+2;
      begin_key = end_val+1;
    }else if (contains("-0123456789.",data[begin_val])){
      end_val = data.find_first_not_of("-0123456789.e",begin_val)-1;// allow negative numbers, decimals and scientific notation
      begin_key = end_val+1;
    }else if (data.substr(begin_val,4)=="true"){
      end_val = begin_val+3;
      begin_key = end_val+1;
    }else if (data.substr(begin_val,5)=="false"){
      end_val = begin_val+4;
      begin_key = end_val+1;
    }else{
      throw; // could not detect the type
    }
    val = data.substr(begin_val,end_val-begin_val+1);

    map.insert(std::pair<std::string, std::string> (key,val));
  }
  return map;
}

template<int size>
Eigen::Matrix<double, size, 1> json_to_eigen_d(const std::string& data){
  auto string_list = parse_json_list(data);
  Eigen::Matrix<double, size, 1> vec(string_list.size());
  for (unsigned int i=0; i< string_list.size(); i++){
    vec[i] = parse_json_double(string_list[i]);
  }
  return vec;
}

// TODO: write parse_json_int
// template<int size>
// Eigen::Matrix<int, size, 1> json_to_eigen_i(const std::string& data){
//   auto string_list = parse_json_list(data);
//   Eigen::Matrix<int, size, 1> vec(string_list.size());
//   for (int i=0; i<string_list.size(); i++){
//     vec[i] = parse_json_int(string_list[i]);
//   }
//   return vec;
// }

double parse_json_double(const std::string &data){
  // TODO add support for NaN and any other special numeric types
  if (data == "Infinity")
    return std::numeric_limits<double>::infinity();
  else if (data == "-Infinity")
    return -std::numeric_limits<double>::infinity();
  else if (data == "NaN")
    return NAN;
  else
    return std::stod(data);
}

std::size_t find_closing_delim(const std::string &data, std::size_t pos, char open_delim, char close_delim){

  auto closing_pos = data.find(close_delim,pos+1);
  auto opening_pos = data.find(open_delim,pos+1);

  while(closing_pos > opening_pos){
    pos = find_closing_delim(data,opening_pos,open_delim,close_delim);
    closing_pos = data.find(close_delim,pos+1);
    opening_pos = data.find(open_delim,pos+1);
  }

  if (closing_pos == std::string::npos){
    // bad JSON format: there is no corresponding closing bracket
    throw(1);
  }

  return closing_pos;
}

bool contains(const std::string &data, const char* a){
  return data.find(a) != std::string::npos;
}

bool contains(const std::string &data, char a){
  return data.find(a) != std::string::npos;
}

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake