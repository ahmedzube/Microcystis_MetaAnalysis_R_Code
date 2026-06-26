##After creating a new working directory
##Use getwd to confirm the directory of the newly created one
getwd()
##Check the files in the new directory, most likely non since its new
dir()
##Create a new folder data(raw and derived data)
Folder<- "data"
dir.create(Folder)
subfolder<- c("raw-data", "derived-data")
Paths<- file.path(Folder, subfolder)
dir.create(Paths[1])
dir.create(Paths[2])
##Confirm the folders created
dir(Folder)

#########################################################
##########If we're saving a link from online#############
#########################################################
# url  <-  " http://ecor.ib.usp.br/lib/exe/fetch.php?media=dados:caixeta.csv" 
#Paths <- file.path(paths[1], "caixeta.csv")
#download.file(url, destfile = Paths)


###Create folder for code
dir.create("Code")

##Create folder for result
dir.create("Result")

dir()
