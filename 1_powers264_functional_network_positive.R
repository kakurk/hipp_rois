#Author: Arnab Roy, Ph.D. Binghamton University, NY-13902
#Post-Doctoral Researcher, Penn. State University, State College-16801
#Date: 24 July 2015

#----------------------------------------------------------------------------------------------------------------------------------

#Objective: The objective of this code is to make connectivity network using powers 264 ROI

#*****************************************************************************************************************************


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

corr_matrix <- function(sig_file,op_name,p_threshold) #corr_matrix.begins
{

library("lsr")

 cormat <- correlate(sig_file,corr.method="pearson",p.adjust.method="fdr")
 cormat$correlation[is.na(cormat$correlation)] <- 0 #for NA cases set correlation to 0
 cormat$p.value[is.na(cormat$p.value)] <- 10000 #For NA cases set the p-value to very large number

 #The NA cases will anyways get excluded as p-value will be large, we will only choose for p-value < 0.05
  
 #-----------------------------------------------------------------------------------

 #Save only the edges that qualify the p-threshold
            
 op_table <- vector('list',(nrow(cormat$correlation)*nrow(cormat$correlation)+1))
 op_table <- list(c('vox_col_id_A','vox_col_id_B','weight','p_value')) #this is the file header

 counter <- 1 #initialize the counter to 1

 for(c4 in 1:(nrow(cormat$correlation)-1)) #For.c4.begins
    {
     for(c5 in (c4+1):ncol(cormat$correlation)) #For.c5.begins
        {
         if(cormat$p.value[c4,c5] < p_threshold)
           {
            vox_col_id_a1 <- c4
            vox_col_id_a2 <- c5
      
            if(cormat$correlation[c4,c5] > 0)
              {
               op_table[[counter+1]] <- c(vox_col_id_a1,vox_col_id_a2,cormat$correlation[c4,c5],cormat$p.value[c4,c5]) 
              }
            counter <- counter + 1

            }

        } #For.c5.ends

    } #For.c4.ends



 #--------------------------------------------------------------------------------------

 file.create(op_name)

 for(filecounter in 1:counter)
    {
     write(c(unlist(op_table[[filecounter]])),file=op_name,append=TRUE,ncol=100)
    }

} #corr_matrix.begins

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


F_powers264_functional_network <- function(IP_f_path,IP_s_path,IP_subject_folder,IP_c_file,IP_bm_file,IP_powers_coordinate,IP_output_file_tag,IP_p_val_threshold)#the main funtion begins
{

#Load the fmri library

library("fmri")

#----------------------------------------------------------------------------------------------------------------

#c1 controls the subjects
for(c1 in 1:length(IP_subject_folder)) #for.c1.begins
   {
    #the functional nifti file
    f_nifti_file_name <- paste(IP_f_path,IP_subject_folder[c1],"/",IP_c_file,sep="")
    print(f_nifti_file_name)
    BOLD <- read.NIFTI(f_nifti_file_name, level = 0.75,setmask=FALSE)
    BOLD_data <- extract.data(BOLD, what = "data")
    time_length <- dim(BOLD_data)[4]
    #-----------------------------------------------------------------------------------------------------------

    #the structural nifti file
    s_nifti_file_name <- paste(IP_s_path,IP_subject_folder[c1],"/",IP_bm_file[c1],sep="")
    print(s_nifti_file_name)
    mask <- read.NIFTI(s_nifti_file_name, level = 0.75,setmask=FALSE)
    mask_data <- extract.data(mask, what = "data")
    threshold_grey <- 0.6*max(unlist(mask_data))
    #-----------------------------------------------------------------------------------------------------------

    #if mask dimensions do not match the data dimension exit

    if(BOLD$dim[1] != mask$dim[1] && BOLD$dim[2]  != mask$dim[2] && BOLD$dim[3]  != mask$dim[3])
      {system("echo \"mask dimensions not equal to data dimensions\" >> error ");quit();}


    print(c(BOLD$dim[1],mask$dim[1],BOLD$dim[2],mask$dim[2],BOLD$dim[3],mask$dim[3]))
       
    #-----------------------------------------------------------------------------------------------------------

    #first read powers coordinate file
 
    p_cood <- read.table(IP_powers_coordinate,header=TRUE)
  
    #define empty array for storing the ROI signals
    signal_data <- vector('list',278)  #PARA

    #now iterate through all 264 powers coordinates and find out which voxels belong to each coordinate.
    for(c2 in 1:278)#for.c2.begins  #PARA
       {

        #find the voxel closest to this power's coordinate
        vox_x_cood <- round((91 - p_cood$x[c2])/3)
        vox_y_cood <- round((107 + p_cood$y[c2])/3) 
        vox_z_cood <- round((89 + p_cood$z[c2])/3) 

        #print(c(p_cood$x[c2],p_cood$y[c2],p_cood$z[c2],paste(p_cood$system[c2])))
        #print(c(p_cood$x[c2],p_cood$y[c2],p_cood$z[c2],vox_x_cood,vox_y_cood,vox_z_cood))

 
        #now check all 27 voxels around it and see which ones belong to grey matter.
        #and then create an averge signal using these voxels

        voxel_count <- 0 
        signal <- rep(0,205)  #PARA
        for(c3x in -2:2)#for.c3x.begins
           {
            for(c3y in -2:2)#for.c3y.begins
               {
                for(c3z in -2:2)#for.c3z.begins
                   {

                     if(vox_x_cood+c3x > 0 && vox_x_cood+c3x <= BOLD$dim[1])
                       {
                     if(vox_y_cood+c3y > 0 && vox_y_cood+c3y <= BOLD$dim[2])
                       {
                     if(vox_z_cood+c3z > 0 && vox_z_cood+c3z <= BOLD$dim[3])
                       {

                        #print(c(vox_x_cood+c3x,vox_y_cood+c3y,vox_z_cood+c3z))
                        if(mask_data[vox_x_cood+c3x,vox_y_cood+c3y,vox_z_cood+c3z,1] > threshold_grey)
                       {
                        voxel_count <- voxel_count + 1
                        signal <- signal + BOLD_data[vox_x_cood+c3x,vox_y_cood+c3y,vox_z_cood+c3z,1:time_length] #PARA
                       }
                      }
                      }
                      }

                   }#for.c3z.ends
               }#for.c3y.ends
           }#for.c3x.ends

        #print(c(c2,voxel_count))

        
        #store the average signal
        if(voxel_count == 0)
          {signal_data[[c2]] <- rep(0,205)}  #PARA
        else
          {signal_data[[c2]] <- signal/voxel_count}

       }#for.c2.ends

    #-----------------------------------------------------------------------------------------------------------

    #now create the correlation matrix
    #please note that the signal for each ROI is store across columns.
    #that is col-1 is time-1, col-2 is time-2, and so on.
    #for correlation to work properly we need to transpose this matrix

    signal_data_transposed <- t(do.call(rbind,signal_data))

    corr_matrix(signal_data_transposed,paste(IP_output_file_tag,"_",IP_subject_folder[c1],sep=""),IP_p_val_threshold)
    print("Correlation matrix created.")



   }#for.c1.ends


#----------------------------------------------------------------------------------------------------------------

}#the main funtion ends

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------



functional_path <- '/mnt/psych-bhampstelab/fhillary/a_data/analysis/5_subsystems/1_createnetwork_positive/hippocampus/for_kyle/conn_files/'

structural_path <- '/mnt/psych-bhampstelab/fhillary/a_data/analysis/5_subsystems/1_createnetwork_positive/hippocampus/for_kyle/resample_files/'

subject_folder <- c('subject_AM05','subject_AM10','subject_AM15','subject_AM22','subject_AM23','subject_AM24','subject_AM35','subject_AM36','subject_AM53','subject_AM54','subject_AM55','subject_AM58','subject_AM60',
'subject_AM64','subject_AM77','subject_AM79','subject_AM80','subject_AM82','subject_AM84')



conn_output_file<- 'conn_project1_new/results/preprocessing/niftiDATA_Subject001_Condition000.nii'  #PARA

binary_mask_file <- c('resample_wc1cr_co20100910_094940t1mprageSAGNEWs014a1001.nii', #AM05
'resample_wc1cr_co20110401_085207t1mprageSAGNEWs012a1001.nii',#AM10
'resample_wc1cr_co20110824_093643t1mprageSAGNEWs012a1001.nii',#AM15
'resample_wc1cr_co20120127_113849t1mprageSAGNEWs012a1001.nii',#AM22
'resample_wc1cr_co20120112_115928t1mprageSAGNEWs012a1001.nii') #AM23




powers_coordinate <- './powers_MNI_coordinate.dat'  #PARA

output_file_tag <- 'powers_264_conmat_MCIpost_1'  #PARA

p_val_threshold <- 0.05 #PARA


F_powers264_functional_network(functional_path,structural_path,subject_folder,conn_output_file,binary_mask_file,powers_coordinate,output_file_tag,p_val_threshold)






