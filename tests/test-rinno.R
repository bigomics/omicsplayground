###https://ficonsulting.com/filabs/RInno

install.packages("RInno")
require(RInno)
RInno::install_inno()
example_app(app_dir = "app") 
create_app(app_name = "myapp", app_dir = "app")
compile_iss()
