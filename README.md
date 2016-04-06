# ExoJulia <img src="https://github.com/jlustigy/ExoJulia/blob/master/Extras/proto1.png" width="13%" height="13%" align="left" />

**Julia packages for fast exoplanet modeling**

This code is being developed by Astronomy graduate students at the University of Washington for the course *ASTR 598A: Exoplanets*. The course and this project are led by [Professor Eric Agol](http://faculty.washington.edu/agol/index.html).  

## ASTR 598A Collaboration Instructions

### Getting started

You should only need to perform the following once:

1. **Fork this repository** (click button in upper right of this page) to create a working branch on *your* github account
2. **Clone *your* fork** onto your computer:
  
  ```bash
  git clone git@github.com:YOUR_USERNAME/ExoJulia.git
  ```

3. Add *this* repo as an [upstream remote for your fork](https://help.github.com/articles/configuring-a-remote-for-a-fork/):

  ```bash
  git remote add upstream git@github.com:jlustigy/ExoJulia.git
  ```

3. Do HW assignments in the `Homework/` subdirectories (see below)

Now you can code on your machine and commit/push changes to your remote github repo without worrying about messing up or being messed up by other peoples work-in-progress code.  

e.g. if you write a module over the course of a few days, it's nice to push to github a few times over those days to ensure you're backing everything up remotely. With git+github the dogs can hang around our homework, but they can never eat it. *Good dog*.  

### Each week
1. [**Sync your fork** with the upstream branch](https://help.github.com/articles/syncing-a-fork/) to ensure that we are all writing code that uses the same ExoJulia package:    

  ```bash
  git fetch upstream  
  git checkout master  
  git merge upstream/master  
  ```  
  This brings your fork's master branch into sync with the upstream repository, without losing your local changes.
    
2. **Create a unique homework directory**: `Homework\hw#\Name1_Name2` where "#" is the homework number, and "Name1" and "Name2" (and so forth) are the names of the group members. So each homework assignment will have a directory, within which each coding group will have a directory. This will ensure everyone's code is saved for future reference and will (hopefully) minimize merge conflicts that can occur when people have worked on the same files.  
3. **Write code** and work in Jupyter Notebooks within the directory you just created. Make sure to reference the "official" (previously selected) code if you are building on previous work. DO NOT LOAD MODULES FROM OTHER HOMEWORK DIRECTORIES.
4. **Add, commit, push** to your remote branch:
  
  ```bash
  git add -A
  git commit -am "Message"
  git push origin master
  ```
  
5. **Submit pull request** by navigating to your fork on github.com and pressing the green button for "New Pull Request". This will request that your fork be "pulled" onto the Master branch. Since you've done all your work in a compeletly new directory (didn't you?) there won't be any issues merging :) 
6. Each week the fastest code (and then cleanist in the case of ties?) will be selected to become the "offical" code in the ExoJulia package.  

## Useful notes

* [Julia progamming style guide](http://docs.julialang.org/en/release-0.4/manual/style-guide/)
