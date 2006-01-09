# This is Google diff-match-patch r103
# http://google-diff-match-patch.googlecode.com/

add_library (GDiff gdiff.h gdiff.cpp)

set_target_properties(GDiff PROPERTIES FOLDER Testing)