##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message("===============================================================")
message("======================= init.R ================================")
message("===============================================================")

## Parse access logs
access.dirs = c("/var/www/html/logs", "/var/log/apache2","/var/log/apache",
                "/var/log/httpd","/var/log/nginx","../logs")
access.dirs <- access.dirs[dir.exists(access.dirs)]
access.dirs
##ACCESS.LOG <- pgx.parseAccessLogs(access.dirs[], filter.get=NULL)
ACCESS.LOG <- pgx.parseAccessLogs(access.dirs[], filter.get="playground")
names(ACCESS.LOG)
sum(ACCESS.LOG$visitors$count)

##-----------------------------------------------------
## Initialize ORCA server
##-----------------------------------------------------
## see: pgx-module.R
ORCA <- initOrca(launch=TRUE) 
class(ORCA)
if(is.null(ORCA)) {
    warning("##### FATAL:: Could not connect to ORCA server. Please start ORCA. #####")
    stop()
}
