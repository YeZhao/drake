#import <iostream>
#import "parse_json.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

int main(){
  // TODO: test scientific notation in the parser
  std::string nested_braces ="{{{}{{}{}}{}}{{}{}}{{}}}";
  std::string nested_brackets ="[[[[][][]][][]]]";
  std::string test_string = "\"This is a JSON string\", other junk";
  std::string test_map= "{\"1\":\"hello\",\"2\":\"world\",\"test\":4,\"test\":true,\"test2\":false,\"asdg\":[1,2,3,4,5]}";
  std::string test_list= "[1,2,3,true,false,\"test\",{\"1\":1,\"2\":2},[1,2,3,4]]";
  std::string test_numeric_list= "[1,2,3,.1,0.2,-1.3,-Infinity,Infinity,NaN,0.0,-0.0]";
  
  // should be true
  std::cout << contains("asdfg","df") << std::endl;
  // should be false
  std::cout << contains("asdfg","dfa") << std::endl;
  // should be true
  std::cout << contains("asdfg",'d') << std::endl;
  // should be false
  std::cout << contains("asdfg",'!') << std::endl;
  // should be 12
  std::cout << "Closing brace position: " << find_closing_delim(nested_braces,1) << std::endl;
  // should be 22
  std::cout << "Closing brace position: " << find_closing_delim(nested_brackets,2,'[',']') << std::endl;
  // should be 9
  std::cout << "Closing string position: " << find_closing_delim(test_string,0,'\"','\"') << std::endl;
  // TODO: check the contents of the output map
  parse_json_object(test_map);
  // TODO: check the contents of the output list
  parse_json_list(test_list);
  std::cout << "Testing JSON to Eigen\n" << json_to_eigen_d<Eigen::Dynamic>(test_numeric_list) << std::endl;

  return 0;
}

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake

int main() {
  return drake::examples::kuka_iiwa_arm::main();
}