# Feedback

Overall the structure of the app is well designed and the code handles the great complexity very well. Very good use of modules and a clear separation between shiny code vs business logic is certainly making our job easier!

**Don't be alarmed by the large amount of feedback, that's not an indication that the code is bad!** Much of the feedback is centered around defensive programming and maintainability - making the application more robust, less prone to errors, and making it quicker to develop.

## Highest Priority

- The local installation needs to be able to work without errors. There are several packages that are being used which are no longer available on CRAN/BioConductor. If these packages are required for the application to run, it would help to install the specific version of the package that works from the archives.
	- There also seem to be a few bugs in the data creation scripts, because even after all the packages have been installed, none of the build scripts run successfully.
	- Even after installing all packages I initially thought I need to install, when opening certain scripts there are packages which still aren't installed such as `rworldmap`.
- I recommend using the [`{config}`](https://rstudio.github.io/config/) package for any settings. Settings should not be hard-coded in the code. This way you can have one version in production, in GitHub, in Docker, in your own local environment. When working locally and some setting has to be different from production, you won't need to remember to always change the setting before commiting.
  - There are several absolute paths in the project, for example `/tmp/...` and `/home/kwee/...`. These should be set in a config file, because currently they only work on your local machine
	- The OPTIONS file settings can also be added to this, which might make it easier to parse than the current format. But you can also keep the OPTIONS file separately from a local config file.
	- The settings in `global.R` (USER_MODE, DEV, DEBUG, WATERMARK, etc) should also be added to a config file
- There is a lot of dead code, which makes it more difficult to traverse the code and understand what's important and what's not. All of this added bloat is hurting readability. I didn't delete the dead code but I think you should, and since this is a git repo then you'll always be able to safely retrieve it if you want.
  - In `global.R` and initialization scripts there are many instances of having a variable name or an expression on a line, without it doing anything - presumably as a way for a dev to test the value when running the script line-by-line.
  - There are `if(0)` statements - those should use some boolean flag instead because it's not clear if it's meant to ever run or if it's test code.
  - There are functions that are defined and never get used (for example `social_buttons()` and `premium.feature()`) and there are functions that are defined multiple times but with different code (such as `dbg()`).
  - There are also chunks of code that aren't used - for example the javascript code for detecting browser dimensions. Shiny never uses it. 
  - If any of the `if (0) {...}` chunks are meant as testing, then it would be sensible to create a testing suite that can be run whenever a change has been made. This way, you can create more robust tests that can check a variety of scenarios so that when users upload custom data, then you are less likely to see unexpected behaviour. {testthat} is the recommended package for testing. Learn more here https://r-pkgs.org/tests.html
- There is a tendency to use `renderUI()` + `uiOutput()` to create all of the UI from the server dynamically instead of from the UI. It's much quicker for the application to have HTML generated from the UI rather than the server. It also helps have better understanding of where elements appear in the UI and makes the server logic much shorter and easier to read.
	- This would perhaps eliminate the need to include `outputOptions(output, "ui_element", suspendWhenHidden = FALSE)` as the server is no longer creating that element.
	- I have re-factored `TcgaBoard` and split into two files, one for the UI and one for server logic, plus removing any code that wasn't being run. Note that this hasn't been fully re-factored, the radio buttons sent to `plotModule` can be moved to the UI once that module has been re-factored, and the inline style has been left in.
- I see that you have a question about whether to use `local=TRUE` or `FALSE`: I do recommend using `local=TRUE` when sourcing, because you don't want the global environment to be dirtied with all the objects, rather you want everything to be available for the caller. Right now the global env happens to also be the environment of the script that sources files so it wouldn't make a behavioural difference, but if the shiny app is called from within a function (for example if you make this project into a package) then the app's environment is no longer the globalenv, and then it will matter.
  - With that said, when there are so many files being sourced that there is no organization in who has access to what, it becomes a bit of a wild west. Everything is accessible from everywhere and it's impossible to know where a variable comes from and who should use it.
    - A quick easy "solution" to add some organization and clarity is to simply not create many global variables from many different files. Variables created within the same file should be grouped together and should be more explicit about where they come from. For example, in global.R, instead of creating variables [OPG, RDIR, WATERMARK] and more, you can create a single `.globals <- new.env()` and the assign the same variables into `.globals`. This way it's more clear where a variable came from, they're grouped together, and you limit the number of objects in the global scope.
    - A better solution, but it would also require more work, is to use a more modular framework (not to be confused with shiny's module feature). It would allow you to explicitly import and export certain objects from each source file, which will improve long term maintainability. There are two options for this: the {modules} package and the {box} package. I personally recommend {box} since it's a much more widely used package and its author is a well established member of the R community with many great contributions over the years.
- In server.R when initializing all the board modules, all of them take the entire `env` object as a sole parameter. This defeats one of the main purposes of using modules - being explicit about your requirements and limiting the scope of objects you can use/modify. Most of the boards only use one or two values (such as env[["load"]][["inputData"]] or env[["expr"]][["selected_gxmethods"]]) so each module should list the individual parameters it actually needs to use and only take those values
- A significant portion of the startup time is taken up by creating the Orca server (~half a minute for me). It seems the only use of it is to download static plotly plots. I would be worthwhile looking into an alternative such as `plotly::orca` and removing the `initOrca` function.
  - There is some issue between shiny and orca on Windows (at least for me) where it doesn't like the temporary file path, but can be resolved by making a temporary file within the project and copying that file to the `file` argument in downloadHandler.  
- When loading objects in the startup, the largest object is `GSETS` which is a 776MB 63k-long list of vectors of gene names. Reading this huge list takes several seconds, is it absolutely necessary? Maybe consider moving this data into eg. a SQLite database where you can reference the required IDs rather than loading the entire list on init. 

## Project Structure


- With an application as large as this, it is worth splitting out the `app.R` file into separate `ui.R` and `server.R` files. The UI and server logic are two separate pieces of functionality and are easier to deal with in separate files.
	- The creation of constants in the `app.R` and `app-init.R` files can be moved to the `global.R` file. This file automatically gets sourced before anything else, so these variables will always be created.
	- The `global.R` file doesn't need to be sourced explicitly 
	- Any function definitions can be moved to a `utils.R` file.
	- This has been done in my PR
- Since there seems to be a lot of business logic that isn't part of the shiny app (most of the code in the `/R` folder), it would be worthwhile to consider converting this project into an R package. All of the non-shiny code would be exported from the package, and then you gain some of the best practices and tools that packages have. It'll naturally force your code to be more robust. The shiny app can be created inside the `inst/` folder and will only contain code that is explicitly related to the web application, and this way it will be easier to manage the codebase - it'll be clear what parts of the code can be developed and tested independently of the shiny app. I don't think it would be a huge mission to convert the current project to a package because it seems the files are mostly organized in such a way that files are either fully shiny related or not. A file that doesn't use `render*()`, `input$`, `output$`, `observe()`, `reactive()` or UI elements is most likely a good file to move into the package and export its functions.
  - If not converting to a package, then since almost all of the files in `/R` `/boards` and `/modules` are being sourced, you can use a loop to source all the files in these directories automatically
- Any file that is deprecated should be removed from the repository, it is preserved in GitHub's version control if you ever need to use it again.
- The majority of the files in `resources` are not being used, or duplicated from the `www`. Try to only include files that are being used.

## Code Style/Best Practices

- With a large project such as this one, I would recommend using Git branches for any piece of development so that the master branch can remain as stable as possible. This way you don't require adding `if (DEV)` throughout the application, and avoids any breaking changes on the master version that are created whilst developing any new module/functionality. This also goes for any file/folder that isn't included in the GitHub repository such as `filesx`.
- This is a large application, so you need to be careful about accidentally overwriting variables and options. I've seen the shiny upload size change from 999Mb to 500Mb, functions being created multiple times with different code, files being sourced multiple times. It's hard to know which is being called last, and this can cause unexpected behaviour.
- There are many `outputOptions(output, "x", suspendWhenHidden=FALSE)` - I strongly recommend removing those, and only keeping the ones that are absolutely required. Shiny has a built in feature where it tries to be smart enough to not render or waste computation time on anything that isn't currently visible on the page, and by disabling suspendWhenHidden you're removing that feature, which makes performance worse.
- Some functions are also masking base R functions, such as `read.csv2`. Avoid ever doing this; R users are familiar with these functions and they may use this function and not understand why it is behaving differently.
	- This also goes for assigning values to names of known functions such as `c`, `max` and `t`. 
	- For shiny applications, avoid using `output` as a variable name, to not add confusion on whehter this is the shiny output object
  - `F` has been used both as a function argument and as a variable too. This is the shorthand for `FALSE` and certain text editors will highlight it as such. It is one of the reasons to always write the longhand `TRUE` and `FALSE`, but also why you should never use `T` and `F` as variable names. In general it's a good idea to avoid variable names with less than 3 characters. 
- There are many warnings about partial argument matching, which need to be fixed. (1) It will remove these warning messages from the console, and a cleaner console means it's easier to notice true warnings/errors/unexpected messages. (2) It's safer to use the full argument name because if the function changes in the future and there is a new parameter matching the shortened argument, then it will use that argument rather than the original argument.
- Several packages are loaded within functions, using both `library` and `require`. (1) You should never use `require` anywhere, you should only use `library` (see [here for explanation](https://stackoverflow.com/a/51263513/3943160)). (2) Inside a function you should use neither, because it makes the function have side effects. Inside functions, use `::` or use `library` in the beginning of the app.
- `%<a-%` has frequently been used to assign reactive variables. The way reactive variables work is that whenever an input or reactive variable used within the expression is changed, the whole expression is run again. This is essentially duplicating the active binding. Using `<-` should suffice. There shouldn't be a need to use active binding for reactive variables.
- Arguments for functions (including modules) should be as simple to use as possible. Sending named vectors or lists as arguments means that they have to be designed in a specific way, otherwise the function will break. And it means that the receiver needs to know the structure of that parameter and its names. By splitting them out into separate arguments, it is easier to define the arguments and know what format the inputs should be.
	- One example of this is the `env` object that is sent as an input to every module. If the `env` object ever changes, currently all of the modules will need checking to see if anything breaks. By being more specific, this can all be adjusted in the top level server script. Currently most of the modules only seem to require `env$load$inputData` and therefore that should be the object sent to each module.
	- Another is the `max.limit` for LoadingBoard (l#32). By splitting out into 4 separate arguments, it is easy to change just one of these without having to create a new vector with all 4 limits.
- When writing functions for use in shiny applications, avoid using the reactive objects as arguments, and instead use the returned value from the reactive object. When using a reactive object in a function, you have to use `()` whenever it is called, and outside of the shiny context, it looks like a function call. This prevents the function from being used outside of the application.
- Try to make variable names a more descriptive. The longer a function is, the more scrolling up and down it takes to find where a variable is defined, and it is a lot easier to know that `dir_name` is a directory path compared to `dd`.
- Having strict coding standards will help to no end with the readability of the application code. The Tidyverse style is a good start, and suggest having a read through it (https://style.tidyverse.org/). Some styling points to focus on:
	- Variable names - try and keep a consistent naming convention. Try to avoid including `.` in variable names; although base R is full of them, the main use of the period in variable names is to create S3 methods for functions, so the function can easily handle different input classes.
	- `{}` - `{` should be the last character on the line. `}` should be the first character on the line. 
	- Assignment - Stick to one method of assignment, ideally `<-` rather than `=` due to the latter also being used in function calls.
	- Quotations - Similar to assignment, either stick to `"` or `'`. `"` is generally regarded as the better option; `'` is used a lot more in text.
	- `;` - Try to avoid using `;` when writing a series of short statements. It might reduce the number of lines used, but it is easy to miss when reading through code. The legibility of code is much more important than saving a couple of lines in a script.
	- Comments - Only one `#` is required for a comment, there is no advantage to adding any more prior to a comment.
- Adding prefixes such as `pgx.` to a multitude of functions reduces the meaningfulness of the prefix. Using the naming conventions mentioned above, these should either be included after an initial verb if relevant to the function, or deleted altogether.
- There are a lot of `for` loops used throughout the application where each iteration creates a value or list that is assigned to an element in a previously empty object. A more efficient way of doing this is by using the `apply` family of functions. These will apply a function to every item in a given array/list. Two in particularly useful are `vapply`, which will return a vector of a specified length and class, whilst `lapply` will return a list. 
	- I have created an example adapting how the navbar panels were being created in the UI. You can see how you can avoid all of the extra assignment of variables with this method. Have a read of https://nicercode.github.io/guides/repeating-things/ for more examples of the `apply` family.  
- Avoid using `ifelse(x, y, z)` when the input is a single value. Only use `ifelse()` when the condition is a vector and you expect a vector answer. Instead, use `if(x) y else z`. There's actually a very big difference between the two, and it isn't always visible. Because `ifelse()` expects a vector input, it returns a vector of the length of its input. This means that if the input is a single value, `ifelse()` will only return a length-one vector. This is why `ifelse(TRUE, 1:5, 10:15)` only return `1` instead of returning `c(1,2,3,4,5)`. The correct version is `if(TRUE) 1:5 else 10:15`.
- There is a fair amount of code duplication. Any time you have to copy and paste more than 3 lines of code, then you are better off creating a new function and using it in both places. Having it in a single function means less code to maintain and update, and it is easier to test the small function.
	- Several modules have small chunks of a reactive expression where 2/3 variables are being applied to the same functionality. Instead of copying and replacing with the new variable, also create a function for these so that you only have to copy the function call line.
- For functions where a default argument is `NULL` and later on if it is the default `NULL` then it is assigned another argument, you can change the default for that argument to the default value. For example:
  ```r
  example_module <- function(input, output, session, param1, param2 = NULL) {
    if (is.null(param2)) param2 <- param
  }
  ```
  can be written as
  ```r
  example_module <- function(input, output, session, param1, param2 = param) {
  }
  ```
- Try to avoid using internal functions in the application. These are designed to not be seen, and aren't always going to be stable. It is safer to create your own function with the required piece of code so that if the internal function changes, you still have a working version.
	- However, for `plotly:::to_JSON("")` just remove it. Normally I would recommend `jsonlite::toJSON` but in this case an empty body is not adding any information to the HTTP request and therefore not required.
- When running HTTP requests, there should be a check for the HTTP status, as the expected behaviour will not happen if the HTTP responses wasn't successful. Use `httr::status_code` or `httr::stop_for_status` to add error handling when the response is not a success (any status code 300 or higher). Any HTTP call can fail for a number of reasons and that must always be checked.

## Code Deletion

The codebase for this application is very large, and although the majority of it is required for the application to run, there is a substantial amount that can be deleted without affecting the application, either because it is never called or it simply does nothing to change the environment. The GitHub repository has a good level of version control, so if any of the deleted code is required in the future, it can be searched and re-added. By removing the points mentioned below, the code will become a lot clearer and easier for users to read and diagnose any issues that may exist.

- You have written a lot of comments, and whilst this is better than no comments, comments should generally be used for explaining why a certain piece of logic is used rather than what the code is doing. With clear function and variable names, someone new to the code will be able to understand what is going on without needing comments, whereas they will only get the why from well written comments.
	- Comments that are questions should be written as GitHub issues, that way they get more visibility. 
	- Some of the comments can also be quite misleading. For example, in the BiomarkerBoard, there is a comment of `## top 100` referencing the top `10 * NFEATURES` rows. Looking at the rest of the expression, this actually comes out to 600.
- I see there a lot of informational/debugging `message()` statements. All of these should be using the `dbg()` function so that you can turn them off, to allow you to see a clearer console sometimes to notice actual issues
	- There are many messages that are not informative and should be removed to not clog up the log files
- Any line which is simply printing a variable (with or without `print()`) can be removed. This is useful when running a line-by-line debugging, but it does nothing when a function is run.
- Variables that are created and never used can also be removed. Some of them might have been in use before code got commented out or deleted, but it just adds extra code to be maintained that isn't required.
- There are also variables that are defined, and then redefined on the next line. This is usually when a shiny input is being used. If you want to test a reactive expression, create a function that runs the business logic, and create tests to check it works as expected. Then use that function for the input and you can delete the original assignment.
- I have found there are several functions that are either labelled as old/moved/notworking or never used in the application (or build process). These should be removed.
	- On a similar note, there is a fair amount of commented out code.
- Any `if` condition that contains `0` or `FALSE` should be removed entirely, and any `if(1)` or `if(TRUE)` should be replaced with just the code inside, without having an if statement

## UI

- The vast majority of customised styling is included in .css files, however there are a few instances where `tags$style` or the `style` parameter is being used - all CSS should be provided in a CSS file.
- Similarly, don't use tags$script with hard-coded code inside - use a JavaScript file.
	- In the case of the dimension code, I don't see `input$dimension` called within the application and so can be removed entirely. 
- When adding styling to elements in a module, try to avoid using the element ID to set customised style, and instead create a new class and use that for the custom styling. If the module ID ever changes, the styling linked with the ID will no longer work, whereas the class will still show as expected.
- There are 16 instances of the ID `sidebar`, one for each tab page. I understand that there is only one sidebar visible at any given time, but they all exist on the page (just not visible) at the same time, which is not legal. You cannot have more than one element with the same ID. 
	- If you pass the namespace to `tabView` then you can set unique IDs for each sidebar.
- Another duplicate element is `data_phenoclustsamples` in the DataViewBoard. One of these needs to be changed to preserve the uniqueness of inputs.
- Any JS script or stylesheet that is added to the UI should be inside `tags$head`.
- Instead of writing HTML as a string and using the `HTML()` function, consider using the `tags` list. This list contains every possible HTML tag including `a`, `b` and `em`, and can avoid a lot of paste calls and long lines.

## Server 

- Avoid adding `output$<id>` render functions in any observe or reactive call. These should all be at the top level of the server/module script, otherwise the output will not be bound correctly and can cause unnecessary errors. In the case of the authentication output situation, it is wrapped by an `observeEvent`, therefore it can be changed to a `eventReactive` and assign `output$login_warning` the returned value of that.
- `reactive` calls within shiny are lazy; that means the expression within the call will only be executed if the reactive object is required elsewhere within the module/server. If an output isn't visible in the UI, then the output rendering will never be called and neither will its dependencies. Because of this, there isn't a need to add conditionals checking whether a reactive should be created or not, as if it isn't required then it won't ever be executed.
- Calling reactive objects cannot be done outside a reactive expression i.e. directly in the server call. There are a few situations where this is done within `if` sections, so it likelt is never run, but this will break the application if run. 
- There are some non-reactive variables that are being created in the server, such as `max.limits`. These can be created outside of the server, as the global environment can be shared across all sessions, meaning it will reduce the load up time for subsequent connections.
- Using `isolate` within an `eventReactive` is unneeded, since the only time the expression will be executed is when the trigger is called. 
- When using `req(x)`, you don't need to include `if (is.null(x)) return(NULL)`. What `req` will do is check whether the elements included are truthy, and if they aren't stop the running the expression and prevent any other reactive/output that uses the reactive the `req` is in to run.
	- To take this one step further, we can use `validate(need(!is.null(x), "message if condition is false"))` to provide a more verbose reason to the user as to why a plot hasn't been generated.
- There are reactive expressions in the server that are several hundred lines long. These could greatly benefit from being transformed into functions outside of the server code. This way they can be chunked up into smaller functions and tested outside of the shiny application.
  - This also goes for any function that is already defined within the server-side code. When in shiny modules that can only be tested manually, whereas when kept in a separate script, they can be used and tested outside of the shiny application.
- Most of the application `reactiveVal` is used when updating counts, however there a couple of time where  `<<-` is used. The issue with this assignment is that it happens in the global environment, and if this is in a module, if the module were reused, this would affect the value in all modules, causing some unexpected behaviour. If this is run on a server, then it will also change the value in all instances. Never use `<<-` (unless you really really know what you're doing and know that you need to).
- Whilst the namespacing in modules is generally well observed, there are a few cases where the id does not include `ns()`.
- Any file that is uploaded by the user should be validated individually, and feedback given to the user about any issues of the uploaded data. Although there is a method for validating the 3 csvs, it is all within the same reactive expression and the return of the validation is in the same expression as the plot generation. By splitting it out into 3 separate `eventReactive` expressions, you can read and validate each file separately and provide clear criteria for each file to pass, then send through a `validate` call to separate it from the plot generation.

## Modules

- Try and keep objects that are only required for one module within that module. For instance, the Firebase authentication objects are created in the server script, but are only used in the authentication module. You can make the FirebaseAuthenticationModule more isolated by initialising within the module, and removes the dependency on the server script. 
- Nearly every page of the application contains a `tabsetPanel` of several pages, all of which are contained in a single module. To reduce the module size into more manageable scripts, you can create sub-modules, where each module contains one of the tabs.
- In modules, try to avoid using `stop` when validating conditions in `observe` statements. This will shut down the whole application when only one module isn't working correctly. The user will not know why the application has crashed without looking at the logs. Instead add some error handling, which may involve a message to the user saying a potential input isn't correct.
- In ClusteringBoard, there is a comment saying the heatmap does not like single gene. Rather than creating a hack around this, create some sort of validation and a clear message to the user that single genes won't produce a heatmap. You can skip a large amount of functionality from this and it is nicer for the user than seeing potentially inaccurate data.
- The `plotModule` is impressive for the amount of plot types that can be passed to it, however with the amount of arguments and all of the `if` clauses to switch between each plot type, it will only work for plots that are currently generated for this application. To make it more generic and accessible for future plots, anything that requires `if` clauses, such as the rendering functions, should be stored as parameters. That way, if a new plot type is introduced, you only have to change the call to the module rather than having to adapt the `plotModule` for the new plot type.
	- By sending the output and render function as function parameters (instead of strings), you avoid the need to write `eval(parse(text = outputFunc))` and can just use `outputFunc`.
	- Sending through the reactive plot as `func` is a bit mis-informational. It should be given a name more accurate such as `reactive_plot`.
- Some of the other arguments for `plotModule` and `tableModule` are rather confusing with no documentation. For example, there are two heights and widths given, and only by looking further into the code we can see why two are required. It is more intuitive to have two height arguments, `plot_height` and `popup_plot_height`. This adds clarification as to why two heights are required and can include the default heights in the function arguments rather than checking the length of `height` in the module.
- Is there a reason in the plotModule to store the file names within the module and not create during the file download process? The parameter `file` in the `content` argument contains the temporary file that can be used to add the exported file to, and changes each time it is called.

## UX

- A lot of the initial loading time is spent loading packages. By using namespaces (`pkg::function` notation) instead of reading the whole library you will significantly reduce loading time (ie. if 10 packages take 2 seconds each to load and they're only used in one specific place in the app, then instead of spending 20s at initialization, you'll spend 2s later on but only if it's required). It will also remove the large list of function masking that occurs at the start of the application.
- While it's nice to welcome users to the application, having two separate modals every time you open the application doesn't provide a great user experience. They have to click twice before even starting to use the application. Feel free to keep the initial modal, but there is no added benefit in my opinion to the second modal. 
- Both the `b` and `strong` tags are being used throughout the application. When looking at the UI of the application they appear the same. However from an accessibility point of view, they slightly differ. Visual readers will read text in `strong` tags in a louder voice to emphasize what is written. Therefore `b` should be used for style, and `strong` should be used for context. This difference is the same for `i` and `em` tags.
- There are a lot of colour schemes used throughout the application. It might make a nicer experience to reduce the number of palettes down to create consistent colour theming across the pages. These can be stored in a separate script so that they can be referenced in any module.
- When uploading a user's own datasets, it only mentions the csv option, it does not mention that pgx files can also be uploaded. It also requires the user to upload all 3 files at the same time, and the information about each file is in a separate box. If you write all of the information as part of the description then it will be clearer to the user what is required.

## Minor Points 

- When writing `for` loops, you don't need to create the iterator beforehand, the loop will assign the first value to the iterator.
- Rather than using `[.]` for special characters in regular expressions, use `\\.` instead. 
- When there is only one element in `c()`, you can avoid writing `c()` and just write the variable.
- I have seen some boards where the height in the UI part of the module is different to the height in the server part of the module. These should be consistent with each other.
- You have a lot of custom `tipify` custom functions (including 2 identical functions), and then use the standard `tipify` as well as the custom versions. Instead, I would add the option list as an argument with a default parameter of the `list(container = "body")` so that if it isn't needed you can replace it, and solely use the custom functions.
- Sometimes `is.nan(X)` is used right after `is.na(X)`. `is.na(NaN)` will return a true value, therefore the second check of NaN is not required.
- Personally for strings that are longer than 120 characters, I like to split down into several lines and paste together. That way the whole string is available on the screen without having to scroll horizontally.
- Avoid referencing data.frame columns by integer where possible, and use the column names instead. The column names are a lot less likely to change compared to column position, and therefore less code maintenance is required. It's also clear for anyone reading the code which column is being referenced.
- This is a rather pedantic point, but the convention used in the sidebar labels is different to that of the dropdown labels. One has capitalisation and the other doesn't, and would look better if it is consistent.
- Rather than using `x == TRUE`, use `isTRUE(x)` (similarly use `isFALSE(x)` over `x == FALSE`). This is because it is a more thorough check and will always provide `TRUE` or `FALSE`. `NA == TRUE` will return `NA` whereas `isTRUE(NA)` will return `FALSE`. 
- There is a slight difference between the plotModule buttons and tableModule buttons, where the tooltip doesn't appear in the tableModule zoom button but does in the plotModule.
- Some of the numbers in the clustered heatmap are outside of the plot panel, and therefore it is unknown what these numbers relate to.

## Specific Points

- Currently the javascript code for determining browser dimensions isn't being used, but if you do want to get the browser dimensions in the future, instead of that code you can use the `{shinybrowser}` package
- The `Sys.setlocale` calls in `global.R` do not work in Windows, so you can check the OS before running
- I can see `useShinyjs()` is called in the UI, but it is also in `pgx-include.R`, and in a non-shiny context it doesn't do anything.
- There is a double sourcing of `pgx-functions.R` and `pgx-files.R` - in `global.R` and `pgx-include.R`
- In `pgx-files.R`, use `trimws()` as a more efficient way of removing outer whitespace
- `pgx.initDatasetFolder` in `global.R` uses verbose = 1, whereas the functions inside use `TRUE`. Aim to be consistent in the variable classes.
- The `header` parameter of `navbarPage()` is meant to accept UI that will be shared among all the tabs. It's not meant to take arbitrary "setup" HTML. Instead, the entire `navbarPage()` can be wrapped in a `tagList` along with the other setup HTML. This also eliminates the double-tagList of the `header` and `footer` variables.
- For the `BiomarkerUI`, there is only one panel in the `tabsetPanel`, and then wrapped by a `fillCol`. When there is only one child, you are adding redundant HTML tags, making the code more complicated than it needs to be.
- CorrelationBoard (l#1229-1232), rather than using several `ifelse` statements to obtain the `cex` value, you can use `findInterval`. With this method, it is a lot easier to add in extra cuts. 
```r
cex_levels <- c(1.2, 0.8, 0.5, 0.2)
dim_cuts <- c(0, 40, 100, 200, Inf)
cex <- cex_levels[findInterval(ndim, dim_cuts)]
```
- DataViewBoard (l#260), a more efficient way is to use `c(head(order(-rho),15), tail(order(-rho),15))`. This way you don't need to reorder as it is already in the required order.
- DataViewBoard (l#647), With the exception of `NULL`, if originally `eg` is of length 0, then when choosing the first element of `eg` will give you `NA`. Therefore only the null check is required.
- DataViewBoard (l#1503), Rather than adding extra empty columns, it might make it look neater to change the width of the table.
- AuthenticationModule (l#1186), you can store the messages in a csv and choose a row at random. It is easier to add to a csv file than updating a module.
- UploadModule (l#647), if `!is.null(uploaded[["pgx"]])` then by default `pgx` is in the uploaded object, therefore the first check is obsolete.
- pgx-modules (l#17), include `ignore.stderr` as well as `ignore.stdout` as it will produce an error message on Windows machines. Also `return` should have brackets afterwards.
- `in.shinyproxy` should be a constant rather than a function, the environment variables aren't changed during the application running, therefore this value won't change.
- Rather than loading, saving and copying the pgx file (LoadingBoard l#240), consider copying it, that should do the same thing and save a lot of time.
- In LoadingBoard, the reactive value `currentPGX` is only assigned to once in an `observeEvent` and then assigned in a `reactive` expression. With the exception of a modal being shown/hidden, there is little else being affected in the server, so should ideally be written as an `eventReactive`, making the code less complicated.
- LoadingBoard (l#593), `unique` on a new character vector, which seems redundant when you are creating the vector in the same call. 
