//sequentially opens all the files in source directory,
//computes the max projection
//then saves the result in destination directory
//under the same name with the suffix "-max" added


macro "batch convert nd2 file to tif format"{
	showMessage("Choose Source Directory");
    dir1 = getDirectory("Choose Source Directory ");
    showMessage("Choose Destination Directory");
    dir2 = getDirectory("Choose Destination Directory ");

    name_inc = getString("File Name Contains", "");

    list = getFileList(dir1);
    setBatchMode(true);
    for (i=0; i<list.length; i++)
    {
        showProgress(i+1, list.length);
	if( (indexOf(list[i],".nd2") != -1) || (indexOf(list[i],".ND2") != -1))
	{
             	if( indexOf(list[i],name_inc) != -1)
		{
			//open(dir1+list[i]);
		    run("Bio-Formats Windowless Importer", "open=["+dir1+list[i]+"]");  
              
            n = indexOf(list[i],".");
			strname = substring(list[i],0,n);  // get the name of the file
			  
            saveAs("tiff",dir2+strname+".tif");	//save as tiff 
       		close();
		}
	}
    }
    setBatchMode(false);
}


