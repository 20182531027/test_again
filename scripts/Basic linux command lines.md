1. Directory operation commands (pwd, cd, mkdir, ls )
2. File reading, creating and editing (less, cat, head; touch; vi/vim)
3. Common scripts operating lines

## Directory operation commands 
### `Where am I?`
```bash
pwd
# 'pwd' stands for 'print working directory'.
```
![[Pasted image 20241020163639.png]]
 
###  `What's in the directory?`
ls which is short for listing. This command will list the contents of the current directory
```bash
ls [options] [directory_or_path]
```

```bash
ls 
# list the file or directory concisely
ls /home/chenqi5
```
![[Pasted image 20241020164037.png]]

```bash
ls -l
# Besides the name of the file/directory, the owner, the creation time and the permissions are shown.
ls -l /home/chenqi5
```
![[Pasted image 20241020164113.png]]


![](https://raw.githubusercontent.com/sudheerj/Linux-cheat-sheet/master/images/permissions.png)

```bash
chmod 755 chenqi5.txt
# usually use 755, all permissions for owner of the file, while others could only read and executable. if you don't want others see your file, you could choose 711
```

```bash
ls -l /home/chenqi5/test
```

`First column`: The first column displays the file type and permissions. Permissions are represented by 10 characters.

- The first character indicates the file type:
    - `-` represents a regular file.
    - `d` represents a directory.
    - `l` represents a link file, typically a symbolic link.
    - Other characters such as `b`, `c`, `s`, `p` represent block device files, character device files, sockets, and FIFOs (pipes), respectively.
- The following 9 characters represent file permissions, divided into three groups, each with 3 characters:
    - The first group (from left to right) indicates the owner's permissions.
    - The second group indicates the group's permissions.
    - The third group indicates permissions for others.
    - Each character can be `r` (readable), `w` (writable), or `x` (executable), or `-` to indicate the absence of that permission.

`Third column`: The third column displays the username of the file's owner.

`Fourth column`: The fourth column displays the group name to which the file belongs.

`Fifth column`: The fifth column shows the size of the file in bytes. For directories, this number is usually 4096 bytes, which is the space occupied by the directory entries, not the total size of the files within the directory.

`Sixth column`: The sixth column displays the date and time when the file was last modified.

`Last column`: The last column shows the name of the file or directory.


```bash
ls -al
# On the basis of ls -l, the hidden files or directories also shown, which beginning with .
```

![[Pasted image 20241020164910.png]]

Feel free to go through the different options and play with them
```bash
ls --help
```

### `How to enter other directory?`
```
cd [path_or_directory]
# cd means 'change directory'
```

```bash
cd ../
# enter the upper folder

cd ~
# enter your home directory
```

`Absolute vs Relative Paths`
two ways to acess the file in chenqi5 in your home file 
```
cd ../chenqi5/file 
# the Relative Paths
```

```
cd /home/chenqi5/file
# the Absolute path, recommend this path
```

### `How to create a new folder?`

```
mkdir [options] directory_name1 directory_name2
```

```
mkdir WGS_analysis
```

![[Pasted image 20241020170657.png]]


```
mkdir WGS_analysis/sample1 WGS_analysis/sample2
# Create two subfolder sample1, sample2 in WGS_analysis folder
```
![[Pasted image 20241020170911.png]]

![[Pasted image 20241020171137.png]]

```
mkdir -p WGS_analysis/sample1 WGS_analysis/sample2
# one step to complete the steps above
# -p, --parents     no error if existing, make parent directories as needed
```

```
mkdir --help 
# to see more options for mkdir
```

### `How to delete folders?`

```
rm -R WGS_analysis/sample2
# the subfolder sample2 in WGS_analysis is deleted
```
![[Pasted image 20241020172242.png]]

```
rm -R WGS_analysis
# the WGS_analysis is removed in current folder.
```
![[Pasted image 20241020172348.png]]

**`please be cautious when using rm`**
## File reading, creating and editing

### `File reading`
```
less chenqi.txt
# input q to exit the reading 
```

![[Pasted image 20241020173208.png]]

```
cat chenqi.txt
# the contents will shown in the terminal directly
```

![[Pasted image 20241020173053.png]]

```
head -1 chenqi.txt
# look the first line of the file 'chenqi.txt'
```
![[Pasted image 20241020175428.png]]

### `File creating`
```
touch self_introduction.txt
```
![[Pasted image 20241020173443.png]]

`It's better to keep no space in the file name, use "_" as linker`
### `File editing`
```bash
vi self_introduction.txt
# click i in keybord,then you can input what you want; if you finish, you can click Esc in keybord, then :wq! to save the modification, or :q! to discard the modification, and last click the Enter to exit.
```

![[Pasted image 20241020174137.png]]

```bash
echo "Qi Chen" > self_introduction.txt
echo "badminton,table tennis,eat and sleep" >> self_introduction.txt
```

```bash
mv self_introduction.txt chenqi5.txt
# revise the file name
```
 
![[Pasted image 20241020174650.png]]


### `File moving or deleting`
besides revise the name of file or directory, mv can also move the file to other directory or path
```bash
mv chenqi5.txt file/test_again/
# move the current file 'chenqi5.txt' to the directory file/test_again/
```
![[Pasted image 20241020174958.png]]

![[Pasted image 20241020175130.png]]

```
rm chenqi.txt
```
![[Pasted image 20241020175702.png]]
**`please be cautious when using rm`**

### `File copy`

```
cp source_file /path/to/target_file
```

```
cp file/test_again/chenqi5.txt ./
# copy the 'chenqi5.txt' in the file/test_again folder to the current folder
```
![[Pasted image 20241020180017.png]]


## common scripts operating lines
```
Command script [input] [output]
```
`for example:`
bash test.sh -i [file/string] -o [file/string] 
python test.py -i [file/string] -o [file/string] 
......
 `In order to avoid bugs, you'd better write the absolute path of scripts, input file and output file.`