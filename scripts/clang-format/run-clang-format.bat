for /f %%f in ('dir /s/b ..\..\*.h ..\..\*.cpp') do clang-format -i -style=file %%f
