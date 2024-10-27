1. RStudio Interface
2. Common data Structues in R
	- 2.1 Vector
	- 2.2 Data frame
 	- 2.3 List
3. Package usage, such as [ggplot2](https://rstudio.github.io/cheatsheets/data-visualization.pdf), [dplyr](https://rstudio.github.io/cheatsheets/data-transformation.pdf), [ggtree](https://www.jianshu.com/p/b37836d71a70) .etc
# [RStudio Interface](https://rstudio.github.io/cheatsheets/rstudio-ide.pdf)
RStudio is an integrated development environment (IDE) that provides a comprehensive set of tools to enhance your coding experience, making it more productive and enjoyable. RStudio is currently one of the most popular IDEs for R.
![](attachments/Pasted%20image%2020241027103947.png)



> [!Tip] 

 Common usage:
 - "Tab" completion to finish fuction names, file paths and .etc
 - help: ?var or search var directly in help interface
 - check the varibles, dataframe, function .etc in Enviroment 
 - find the commands  you have run in History, or by type the  up button in the keyboard

# Get Stared

```r
getwd()
# It's a good habit to know the work directory before coding.

```

# Common data Structues in R
## Vector
A **vector** is the basic data structure in R. It is composed by a series of values that can be numbers, characters, or logical values.
- `c()` is the function to create a vector.
- `<-` is the assignment operator.
- `Variable` is where vector stored.
	 - names must start with a letter, **cannot** start with a number
	 - names **cannot** contain spaces
	 - avoid special characters, such as @, #, etc. '_' is acceptable, for example, my_data.
	 - case sensitive. my_data and My_data are two different variables.
	 - avoid reserved words or existing function names, such as while, if, else, TRUE, etc.

assign a value 4 to a variable `my_var` with the command
```R
my_var <- c(4,5,6,7,8,9,1,2) # define the varible
class(my_var) # print the type of the varible
#  "numeric"
length(my_var) # print the length of the varible

# numerical operation
test <- 1+1
test <- 5-2
test <- 3*4
test <- 8/2
test <- 9 %% 2 
test <- 2^2
test <- sqrt(4)
mean(my_var)
var(my_var)

my_fruit <- c("apple","banana","pear","watermelon","grape","pineapple","cucumber","peach" )

class(my_fruit)
#  "character"

my_like <- c(TRUE,FALSE,FALSE,T,F,T,T,T)
# T and F are abbreviations for TRUE and FALSE, respectively.
class(my_like)
# "logical"
```

suset the vector using `[]`
```R
# subset by index, subset the first three elements in each vectors.
my_var[c(1,2,3)]
my_fruit[c(1,2,3)]
my_like[c(1,2,3)]

# subset with conditions
my_var[my_var >=4]
length(my_var[my_var >=4])

my_like[my_like=="TRUE"] # 
# TRUE TRUE TRUE TRUE TRUE

my_fruit[grepl('^p',my_fruit)] # find the characters begin with p, `^` 
#"pear"       "pineapple" "peach"


?grepl # if you are stranger to any fuction, you can type '?' ahead.

```

## Data frame
### create data frame in R
You can create a data frame using the `data.frame()` function by passing vectors as arguments.
```R
df <- data.frame(fruit=my_fruit,likeness=my_like,number=my_var) # because of the same length of each vector, you can create data frame directly.
```
you can also read the data frame that created in your 
```R
library(readxl)
genes_14 <- read_excel("M_metascape_result.tya0gfkgg.xlsx", sheet = 2)

df1 <- read.csv("all_sample_spatial_plot_info.csv",row.name=1)
```

### see the charateristics of the data frame `df`
```r
# the dimension of data frame
dim(df) # it has 8 in row, 3 in column
# 8 3 
nrow(df)
ncol(df)
colnames(df) # the colnames of the data frame

head(df) # print the first six rows in default
head(df,7) # print the first seven rows 
head(df,3) # print the three seven rows

```

### subset the data frame, also  use `[]`
```R
df[,c(1,2)] # print the first two columns of the data frame

# three methods to print the certain column 
df[,1] # print the first column
df$fruit # friut means the name of the first column
df[["fruit"]]

# print the value by specifying the index of a row or column
df[1,2] # print the value in row 1 column 2
df[1,"likeness"] # print the value in row 1 as well as the "likeness" column

# subset the number bigger than 4
df[df$number > 4,]

# subset the likeness is "TRUE"
df[df$likeness == "TRUE",]

```

### add new row or column:
```R
# add the ninth row
df[9,] <- c("oranges","TRUE","9")
df <- rbind(df, c("oranges","TRUE","9"))

# add the new column "district"
df$district <- c(rep("China",7),rep("Africa",2))

```

### save the data frame
```
write.csv(df,"test_fruit.csv")
```
## List
use `list()` to create list
```
rec <- list(name="Ming Li", age=30, scores=c(85, 76, 90)) # different to data frame list doesn't need the same length
```

```R
# look the elements in the list

rec[1]
rec$name
rec[['name']]

rec$scores
rec$scores[1]

# delete the element
rec <- rec[-2]
```
you could take a [test](https://app.datacamp.com/learn/assessments?technologies=1) to check your learning.
# Package usage, such as ggplot2, dplyr, ggtree .etc
![](https://mmbiz.qpic.cn/mmbiz_jpg/2sn6vnlXuibpcibBxbIXSrv2mfOG3E5QDnUrTq4dKA19R96d0QSicE2n0MeMibu5War80j7RrjRapVXdZ8oznH7vJQ/640?wx_fmt=jpeg&tp=webp&wxfrom=5&wx_lazy=1&wx_co=1)

the mainly useful packages and their roles in data engineering

![](https://mmbiz.qpic.cn/mmbiz_png/2sn6vnlXuibr3exiclzzKDblg76MApDzRicS8VFc6BooxakyYN8fmfxmjVgb7WeZpYKqdrJQhVfIzsWRNtnicOibFWA/640?wx_fmt=png&tp=webp&wxfrom=5&wx_lazy=1&wx_co=1)

- package install
```
install.packages("ggplot2")
```

- [ggplot2](https://rstudio.github.io/cheatsheets/data-visualization.pdf)
```R
library(ggplot2)

# scatter plot
ggplot(data=mpg)+
	geom_point(mapping=aes(x=displ,y=hwy))
```

```R
head(mpg) # mpg is the self-contained data frame
```

```R
ggplot(data=mpg)+
	geom_point(mapping=aes(x=displ,y=hwy,color=class))
```

```R
# line fit
ggplot(data=mpg)+
	geom_smooth(mapping=aes(x=displ,y=hwy))

```

```R
# bar plot
ggplot(data=diamonds)+
	geom_bar(mapping=aes(x=cut,fill=cut))

```

