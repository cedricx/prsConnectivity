
# coding: utf-8

# In[76]:


def geneMNItoSchaefer (donor, numParcels, numNetworks, mm):

    import itertools
    import os
    import subprocess
    import numpy as np
    import nibabel as nib
    from nilearn import plotting
    import pandas
    import pickle
    
    print 'Start to map genes from donor %s to Schaefer at %s parcels and %s networks at %s mm resolution'         % (donor, numParcels, numNetworks, mm) + '\n'
    #files
    wk_dir = '/Users/hxia/Desktop/BBL/' #change as appropriate
    project_path = os.path.join(wk_dir,'data/joy/BBL/projects/prsConnectivity')
    
    schaefer_path = os.path.join(wk_dir,'data/joy/BBL/studies/pnc/template/Schaefer/MNI')
    parcel_filename = os.path.join(schaefer_path, 'Schaefer2018_%sParcels_%sNetworks_order_FSLMNI152_%smm.nii.gz'                                % (numParcels, numNetworks, mm))
    community_filename = os.path.join(schaefer_path,'Schaefer2018_%sParcels_%sNetworks_order.txt'                                      % (numParcels,numNetworks))
    
    # Genes
    ABI_path = os.path.join(wk_dir,project_path,'ABI')
    gene_filename = os.path.join(ABI_path,'normalized_microarray_donor%s' % donor,'SampleAnnot.csv')
    
    # Parcellation
    img = nib.load(parcel_filename)
    img_data = img.get_data()
    plotting.plot_img(img,draw_cross=False,title='Schaefer %s Parcels' % numParcels,colorbar=True)

    # Community Assignment
    community = pandas.read_csv(community_filename,delim_whitespace=True,header=None)
    
    # Genes
    gene = pandas.io.parsers.read_csv(gene_filename,delimiter=',')
    gene_mni_coord = gene.iloc[:,-3:]
    
    # Parse Community file
    community_assignment ={parcel[0]: parcel[1].split('_')[2] 
                       for idx, parcel in community.iterrows()}
    community_assignment[0] = 'NotAssigned'
    
    # Initiate an empty dataframe of nrow = num gene, ncol = 3 for x y z coordinates
    gene_vox_coord = pandas.DataFrame(np.zeros(gene_mni_coord.shape))
    gene_vox_coord.columns = ['img_x','img_y','img_z']
    
    print 'Mapping genes to image coordinates' + '\n'
    # Loop through each gene coordinate in MNI space
    # and then apply FSL std2imgcoord to transform the points
    # to the image voxel coordiates

    for ridx, row in gene_mni_coord.iterrows():
        coord_cmd = "echo %f %f %f|std2imgcoord -img %s -std %s -vox "         % (row[0],row[1],row[2],parcel_filename,parcel_filename) # a shell command to apply FSL std2imgcoord
        vox_coord_str = subprocess.check_output(coord_cmd, shell = True) # get the output from the shell
        vox_coord_str_list = vox_coord_str.split()
        for cidx, coord in enumerate(vox_coord_str_list): 
            vox_coord_float = float(coord)
            gene_vox_coord.iloc[ridx][cidx] = vox_coord_float
    gene_vox_coord_rd = np.round(gene_vox_coord).astype(int)
    
    gene_parcel_assignment ={idx: img_data[tuple(gene)].astype(int)
                       for idx, gene in gene_vox_coord_rd.iterrows()}
    
    gene_community_assignment = {idx : community_assignment[gene_parcel]
                            for idx, gene_parcel in gene_parcel_assignment.iteritems()}
    
    # Concat a master loopup table of gene assignment
    gene_assignment = pandas.DataFrame([gene_parcel_assignment, gene_community_assignment]).T
    gene_assignment.columns = ['d{}'.format(i) for i, col in enumerate(gene_assignment, 1)]
    gene_assignment.columns = ['parcel','community']
    gene_assignment = pandas.concat([gene_mni_coord,gene_vox_coord,gene_assignment],axis=1)
    
    print gene_assignment.head(10)
    print '\n'
    #Save the variables
    variables_filename = os.path.join(project_path,'ABI/gene_mapping','%sdonor_%sParcels_%sNetwork_%smm.pkl'                                       % (donor,numParcels,numNetworks,mm))
    with open(variables_filename, 'w') as f:
        pickle.dump([gene_mni_coord, gene_vox_coord, gene_vox_coord_rd,                     gene_parcel_assignment, gene_community_assignment,gene_assignment], f)
    #write gene_assignment to a csv
    csv_filename = os.path.join(project_path,'ABI/gene_mapping','%sdonor_%sParcels_%sNetwork_%smm.csv'                                 % (donor,numParcels,numNetworks,mm))
    gene_assignment.to_csv(path_or_buf=csv_filename)
    
    print 'Saved variables at ' + variables_filename + '\n'
    print 'Saved csv at ' + csv_filename + '\n'

