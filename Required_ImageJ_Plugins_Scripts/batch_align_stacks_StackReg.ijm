//Align stacks for FRAP analysis
//Wulan Deng, April 2019


macro "batch align stacks with StackReg"{
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
        if( indexOf(list[i],name_inc) != -1)
		{
			open(dir1+list[i]);
			run("Grays");
			setSlice(50); // move the frame to be aligned with, in the case of FRAP it is the frame either before or after bleach. 
			run("StackReg ", "transformation=[Rigid Body]");

			
            n = indexOf(list[i],".");
			strname = substring(list[i],0,n);
			saveAs("tiff",dir2+strname+"_StackReg.tif");	
       		close();
		}
	}
    
    setBatchMode(false);
}

