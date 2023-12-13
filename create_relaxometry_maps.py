import os, glob, shutil, tarfile
from datetime import datetime
from datetime import timezone
import numpy as np
from pydicom import dcmread
import pandas as pd
import json
import boto3
import sys
import argparse
from pathlib import Path
import inspect



def unpack_qalas_from_targz(tar_path, output_path):
    '''Function unpacks dicom files from tar and
    return path to QALAS scan
    
    The function will first unpack the targz
    files under the output_path, then it will
    find all terminal folders under the unpacked
    data directory. For each terminal folder, the
    function will aggregate a list of the folder's
    files and then attempt to load the first one
    using dcmread, and look at dcm field 
    (0x0008, 0x103e) for infromation about the
    scan. All folders that have the sequence 'QALAS'
    or 'MAGIC' within this field (does not need to be
    exact match) will then be returned in a list.
    
    '''
    
    try:    
        tarf = tarfile.open(tar_path, 'r:gz')
        tarf_folder = os.path.join(output_path, tar_path.split('/')[-1].split('.')[0])
        tarf.extractall(path = tarf_folder)
        os.system('chmod ug+rwx -R {}'.format(tarf_folder)) #Change permissions so we dont run
                                                            #into deleting issues later
    
        qalas_folders = []
        supplemental_infos = []
    
        sub_dirs = list(set(find_terminal_folders(tarf_folder)))
        for temp_dir in sub_dirs:
            files = glob.glob(os.path.join(temp_dir, '*'))
            if len(files):
                for num in range(min(5, len(files))):
                    tmp_dcm = dcmread(files[num])
                    if ('QALAS' in tmp_dcm[0x0008, 0x103e]._value.upper()) or ('MAGIC' in tmp_dcm[0x0008, 0x103e]._value.upper()):
                        qalas_folders.append(temp_dir)
                        temp_dict = {}
                        temp_dict['SeriesInstanceUID'] = str(tmp_dcm[0x0020, 0x000e]._value)
                        temp_dict['StudyInstanceUID'] = str(tmp_dcm[0x0020, 0x000d]._value)
                        supplemental_infos.append(temp_dict.copy())
                        break
    except:
        print('Error: Unable to successfully parse through the dicoms of the current archive. Skipping processing for this archive.')
        sys.stdout.write('Error: Unable to successfully parse through the dicoms of the current archive. Skipping processing for: {}\n\n'.format(tar_path))
        sys.stdout.flush()
        return None, None
        
    return qalas_folders, supplemental_infos


def convert_single_tar(tar_path,
                       working_dir_base_path,
                       container_path,
                       layout_path,
                       global_path,
                       dcm2bids_config_path):
    '''Function to unpack qalas from dcm, convert to maps then bids
    
    The function takes as input the path to a dicom folder with targz
    extension, and the base path for a working directory. The function
    then looks for any DICOM series with QALAS/MAGIC in its name,
    and it unzips them under a folder created with a timestamp at
    working_dir_base_path. If there are more than one qalas scans
    found, only the second one will be converted to quantitative
    maps using the symri container. Then dcm2bids will convert the
    resulting dicoms to bids. From there the function returns
    a dictionary with info needed to understand what happened.
    
    '''

    output_info = {'tar_name' : tar_path.split('/')[-1]}
    
    #Create working directory with time stamp to store BIDS output
    date_time = datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S")
    working_dir = os.path.join(working_dir_base_path, 'work_' + date_time)
    tmp_bids_dir = os.path.join(working_dir, 'tmp_bids', 'final')
    if os.path.exists(working_dir):
        raise ValueError('Error: working dir at {} should not already exist when running this script'.format(working_dir))
        pass
    else:
        sys.stdout.write('Making working directory at: {}\n'.format(working_dir))
        sys.stdout.flush()
        os.makedirs(working_dir)
        os.makedirs(tmp_bids_dir)

    #Parse the targz file name to find the subject and session label
    [subject_label, session_label] = targz_to_sub_ses_labels(tar_path)
    
    output_info['working_dir_path'] = working_dir
    output_info['bids_path'] = tmp_bids_dir
    output_info['processing_date'] = date_time
    output_info['subject_label'] = subject_label
    output_info['session_label'] = session_label
    
    #If the subject name can't be correctly parsed then don't process
    if subject_label == False:
        output_info['subject_label'] = 'NON_CONFORMANT_TAR_NAME'
        output_info['session_label'] = 'NON_CONFORMANT_TAR_NAME'
        output_info['converted_series_name'] = str(-1)
        output_info['num_qalas_scans'] = 0
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 0
        shutil.rmtree(output_info['working_dir_path'])
        return output_info

    #Unpack ONLY the qalas folders from targz
    qalas_folders, supplemental_infos = unpack_qalas_from_targz(tar_path, working_dir)
    if type(qalas_folders) == type(None):
        output_info['converted_series_name'] = str(-1)
        output_info['num_qalas_scans'] = -1
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 1
        shutil.rmtree(output_info['working_dir_path'])
        return output_info
        
    if len(qalas_folders) >= 1:
        output_info['converted_series_name'] = []
        for temp_folder in qalas_folders:
            output_info['converted_series_name'].append(temp_folder.split('/')[-1])
        output_info['num_qalas_scans'] = len(qalas_folders)
    else:
        output_info['converted_series_name'] = str(-1)
        output_info['num_qalas_scans'] = 0
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 0
        shutil.rmtree(output_info['working_dir_path'])
        return output_info
    
    
    initial_log_path = os.path.join(working_dir, 'symri_container.log')
    
    #Specify the base commands for qalas and dcm2bids that will be configured for each qalas folder
    qalas_base_command = 'singularity run -B {global_path}:/opt/symri/bin/global.ini -B {layout_path}:/layout_path -B {qalas_folder}:/input -B {dcm_maps_path}:/output {container_path} --batch-mode --dataset /input --output /output --layout /layout_path --force-anatomy "infant brain" &> {log_path}'
    dcm2bids_base_command = '{dcm2bids_executable_path} -d {dcm_maps_path} --participant {subject_label} --session {session_label} -c {dcm2bids_config_path} -o {tmp_bids_dir}'
    #dcm2bids_executable_path = '/home/umii/shared/conda_environments/leex6144/imaging_python_39/bin/dcm2bids'
    dcm2bids_executable_path = 'dcm2bids'

    #Run symri container
    dcm_maps_paths = []
    for run_num, temp_folder in enumerate(qalas_folders):
        dcm_maps_paths.append(os.path.join(working_dir, temp_folder.split('/')[-2] + '_qalas_derived_dcm_maps_run-{}'.format(run_num)))
        dcm_maps_path = os.path.join(working_dir, temp_folder.split('/')[-2] + '_qalas_derived_dcm_maps_run-{}'.format(run_num))
        os.makedirs(os.path.join(working_dir, dcm_maps_paths[-1]))
        symri_container_command = qalas_base_command.format(global_path=global_path, layout_path=layout_path, container_path=container_path, qalas_folder=temp_folder, dcm_maps_path=dcm_maps_path, log_path=initial_log_path)
        os.system(symri_container_command)
        sys.stdout.write('SyMRI Command:\n')
        sys.stdout.write(symri_container_command + '\n\n')
        
    #Run dcm2bids conversion
    output_info['num_niftis_generated'] = 0
    for run_num, temp_folder in enumerate(qalas_folders):
        if len(glob.glob(os.path.join(dcm_maps_paths[run_num], '*'))) == 6:

            tmp_bids_dir = os.path.join(working_dir, 'tmp_bids', str(run_num))
            dcm2bids_command = dcm2bids_base_command.format(dcm2bids_executable_path=dcm2bids_executable_path, dcm_maps_path=dcm_maps_paths[run_num], dcm2bids_config_path=dcm2bids_config_path, tmp_bids_dir=tmp_bids_dir, subject_label=subject_label, session_label=session_label)
            os.system(dcm2bids_command)
            shutil.copyfile(initial_log_path, os.path.join(tmp_bids_dir, 'sub-' + subject_label, 'ses-' + session_label, 'anat', 'sub-{}_ses-{}_acq-QALAS_desc-SymriContainer.log'.format(output_info['subject_label'], output_info['session_label'])))

            output_info['symri_conversion_error'] = 0
            output_info['num_niftis_generated'] += len(glob.glob(os.path.join(tmp_bids_dir, 'sub*','ses*','anat','*.nii.gz')))
        
            #Update the metadata in the json files
            json_files = glob.glob(os.path.join(tmp_bids_dir, 'sub-' + subject_label, 'ses-' + session_label, 'anat', '*json'))
            for temp_json_file in json_files:
                update_bids_json(temp_json_file, supplemental_info=supplemental_infos[run_num])

            tmp_bids_files = glob.glob(os.path.join(tmp_bids_dir, 'sub-*', 'ses-*', 'anat', '*acq-QALAS*'))
            for tmp_bids_file in tmp_bids_files:
                orig_file_name = tmp_bids_file.split('/')[-1]
                new_file_name = orig_file_name.replace('_acq-QALAS_', '_acq_QALAS_run-{}_'.format(run_num + 1))
                new_file_path = os.path.join(output_info['bids_path'], 'sub-' + subject_label, 'ses-' + session_label, 'anat', new_file_name)
                if os.path.exists(os.path.dirname(new_file_path)) == False:
                    os.makedirs(os.path.dirname(new_file_path))
                shutil.copyfile(tmp_bids_file, new_file_path)
                
            sys.stdout.write('dcm2bids Command:\n')
            sys.stdout.write(dcm2bids_command + '\n\n')
            sys.stdout.flush()



        else:
            output_info['num_niftis_generated'] = 0
            output_info['symri_conversion_error'] = 1
            sys.stdout.write('SyMRI Container was unable to convert QALAS dicoms to Synthetic Maps\n')
            sys.stdout.flush()
	
    
    return output_info


def grab_from_s3(processed_tgz_files, batch_size = None, bucket_name = None, extension = 'tar.gz',
                    output_folder = None, check_date = True,
                    different_config_path = False):
    '''Function to grab dicoms from s3 bucket for symri
    
    This function takes as input a list of tgz files that have
    already been processed and therefore should not be downloaded.
    For the bucket of interest, the function gathers all tgz files,
    cross references them against the list of processed files,
    and then downloads new files as necessary. The user can
    set specific fields such as batch_size (the max number of
    files to download), the extension to look for in files,
    and the output_folder name to store the downloaded
    s3 files.
    
    Parameters
    ----------
    
    processed_tgz_files: list of strings
        Files that are already processed/should not be downloaded
    batch_size: int or none
        The maximum number of files to download. Set to None to
        download as many as are present in the bucket
    bucket_name: str
        The name of the bucket without the s3:// prefix
    extension: str
        The extension used to screen file names
    output_folder: str
        The path to the base folder to store downloaded files.
    check_data: boolean, default True
        If True, only subjects whose data has been static for
        at least 1 day will be downloaded
    different_config_path : bool or string
        Either False (meaning default config
        path is used) or a string, pointing
        to an alternative config path
        
        
    Returns
    -------
    
    targz_paths: list
        A list of the files that have been downloaded
        
    '''
    
    #If the output folder doesn't exist, then make it
    if os.path.exists(output_folder) == False:
        os.makedirs(output_folder)
        
        
    #Grab config path
    if different_config_path == False:
        config_path = ''
    else:
        if type(different_config_path) != str:
            raise NameError('Error: different config path should eithe be string or boolean')
        else:
            config_path = different_config_path
            
    #Find info from config file    
    with open(config_path, 'r') as f:
        lines = f.read().splitlines()
        for temp_line in lines:
            if 'access_key' == temp_line[:10]:
                access_key = temp_line.split('=')[-1].strip()
            if 'secret_key' == temp_line[:10]:
                secret_key = temp_line.split('=')[-1].strip()
            if 'host_base' == temp_line[:9]:
                host_base = temp_line.split('=')[-1].strip()
                if 'https' != host_base[:5]:
                    host_base = 'https://' + host_base
        
    #Create s3 client
    client = boto3.client(
        's3',
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        endpoint_url =host_base
    )
    
    # Create a reusable Paginator
    paginator = client.get_paginator('list_objects')

    # Create a PageIterator from the Paginator
    page_iterator = paginator.paginate(Bucket=bucket_name)

    #Iterate through bucket to find potential subjects
    bucket_contents = []
    for page in page_iterator:
        for temp_dict in page['Contents']:
            bucket_contents.append(temp_dict)

    #Find files that have the right extension, haven't been modified in the
    #last day, and aren't included in the processed_tgz_files list.
    #If batch size is set, only gather the number of sessions specified
    
    if batch_size is None:
        batch_size = float('inf')
    
    
    keys = []
    counter = 0
    extension_offset = -1*len(extension)
    for temp_object in bucket_contents:

        #Only consider objects that are more than one day old
        if check_date & ((datetime.now(timezone.utc) - temp_object['LastModified']).days < 1):
            continue
        else:
            key_name = temp_object['Key']
            if key_name[extension_offset:] == extension:
                if key_name not in processed_tgz_files:
                    keys.append(key_name)
                    counter += 1
            if counter == batch_size:
                break


    #Download the files
    targz_paths = [] 
    for temp_key in keys:
        temp_local_path = os.path.join(output_folder, temp_key)
        if not os.path.exists(os.path.dirname(temp_local_path)):
                os.makedirs(os.path.dirname(temp_local_path))
        client.download_file(bucket_name, temp_key, temp_local_path)
        sys.stdout.write('Downloading file: {} from bucket {}\n'.format(temp_key, bucket_name))
        sys.stdout.flush()
        targz_paths.append(temp_local_path)
        
    return targz_paths

def push_to_s3(base_bids_dir, subject_label, bucket_name = None,
                    prefix = None, different_config_path = False):
    
    '''Send BIDS folder to S3 bucket
    
    Parameters
    ----------
    base_bids_dir : string
        The directory where bids data is stored
    subject_label : string
        The subject label
    bucket_name : string
        Bucket name to upload data to
    prefix : string
        Prefix to include in front of subject path
        when uploading to bucket
    different_config_path : bool or string
        Either False (meaning default config
        path is used) or a string, pointing
        to an alternative config path
    
    '''

    sys.stdout.write(base_bids_dir)
    sys.stdout.write(subject_label)
    sys.stdout.write(bucket_name)
    sys.stdout.write(prefix)
    sys.stdout.write(different_config_path)
    sys.stdout.flush()
    

    #Grab config path
    if different_config_path == False:
        config_path = ''
    else:
        if type(different_config_path) != str:
            raise NameError('Error: different config path should eithe be string or boolean')
        else:
            config_path = different_config_path
            
    #Find info from config file    
    with open(config_path, 'r') as f:
        lines = f.read().splitlines()
        for temp_line in lines:
            if 'access_key' == temp_line[:10]:
                access_key = temp_line.split('=')[-1].strip()
            if 'secret_key' == temp_line[:10]:
                secret_key = temp_line.split('=')[-1].strip()
            if 'host_base' == temp_line[:9]:
                host_base = temp_line.split('=')[-1].strip()
                if 'https' != host_base[:5]:
                    host_base = 'https://' + host_base
        
    #Create s3 client
    client = boto3.client(
        's3',
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        endpoint_url =host_base
    )
    
    try:
        sys.stdout.write('Uploading data for sub-{}\n'.format(subject_label))
        sys.stdout.flush()
        os.chdir(base_bids_dir)
        files = glob.glob('sub*/ses*/anat/*')
        for temp_file in files:
            response = client.upload_file(temp_file, bucket_name, os.path.join(prefix, temp_file))
    except:
        sys.stdout.write('Either part or all of uploading Failed for sub-{}, the corresponding session will not be save to the log file so that processing will run again later.\n'.format(subject_label))
        sys.stdout.flush()
        return False

    return True

def targz_to_sub_ses_labels(targz_name):
    '''Extract subject/session ids/labels from the targz file names'''
    
    #We expect a bucket file like s3://hbcd-dicoms-pilot/TIUFL0008_880852_V02.tar.gz
    # and want to extract the subject label 880852 and session label V02
    
    sub_label = False
    ses_label = False
    
    split_name = targz_name.split('/')[-1].split('_')
    if len(split_name) >= 3:
        if (len(split_name[0]) == 9) and (len(split_name[1]) == 6) and (split_name[2][0] == 'V'):
            if len(split_name) == 3:
                sub_label = split_name[1]
                ses_label = split_name[2].split('.')[0]
            else:
                sub_label = split_name[1]
                ses_label = split_name[2]
    
    return sub_label, ses_label

def update_tracking_log(original_tracking_log, output_info, upload_status, base_path = None):
    '''Update SyMRI Processing tracking log
    
    Takes original tracking log path, and adds
    tracking info for new subjects (using info
    from output_info). Only subjects who were
    successfully uploaded to s3 bucket will be
    added to the tracking log. Then a new timestamped
    tracking log will be created under the base_path
    folder.
    
    Because subjects without QALAS scans don't have
    any data to upload, they will automatically be
    added to the tracking log.
    
    Parameters
    ----------
    
    original_tracking_log : string
        Path to the tracking log found before processing
    output_info : list of dictionaries
        Contains the metadata used to update the tracking log
    upload_status : list of booleans
        Whether upload was successful
    base_path : str
        Folder where tracking logs are stored
    
    '''
    
    columns = ['tar_name', 'subject_label', 'session_label', 'processing_date', 'converted_series_name', 'num_qalas_scans', 'num_niftis_generated']
    tracking_log = pd.read_csv(original_tracking_log)
    for i in range(len(output_info)):
        if upload_status[i] == True:
            
            temp_dict = {temp_dict: output_info[i][temp_dict] for temp_dict in columns}
            if len(temp_dict['converted_series_name']) > 1:
                temp_dict['converted_series_name'] = ' '.join(temp_dict['converted_series_name'])
            temp_df = pd.DataFrame(temp_dict, columns = columns, index = [1])
            tracking_log = pd.concat([tracking_log, temp_df], ignore_index = True)
            
    date_time = datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S")
    new_path = os.path.join(base_path, 'tracking_log_{}.csv'.format(date_time))
    
    if len(output_info) > 0:
        tracking_log.to_csv(new_path, index=False)
    
    return

def update_bids_json(json_file, supplemental_info = None):
    '''Function that updates the BIDS json for SyMRI Niftis
    
    Removes fields like EchoTime, RepetitionTime and FlipAngle,
    and adds fields related to the SyMRI protocol.
    
    Params
    ------
    
    json_file : string
        The json file to update
    supplemental_info : None or dict, default None
        Used for passing SeriesInstanceUID and StudyInstanceUID
        to the json
    
    '''
    
    fields_to_remove = ['EchoTime', 'RepetitionTime', 'FlipAngle']
    details = {"ImageDescription" : "This is a synthetic image derived from a QALAS scan distributed by SyMRI.\\nQunatitative T1, T2, and PD values are estimated from the QALAS scan using numerical algorithms provided by SyMRI.",
                         "ReferenceDOIs" : ["https://doi.org/10.1016/j.mri.2019.08.031"]}
    
    with open(json_file, 'r') as f:
        json_contents = json.load(f)

    for temp_field in fields_to_remove:
        del json_contents[temp_field]

    if type(supplemental_info) == dict:
        json_contents.update(supplemental_info)
    #json_contents['SeriesInstanceUID'] = 'TBD'
    #json_contents['StudyInstanceUID'] = 'TBD'
    #What other fields to add - mainly wanted series/study instance

    json_contents['DerivativeDetails'] = details
    json_string = json.dumps(json_contents, indent = 4)
    with open(json_file, 'w') as f:
        f.write(json_string)
        
    return

def find_terminal_folders(path, folders = set()):
    '''Function to find terminal folders in directory tree.
    
    Recursively looks for any folders that are (1) terminal
    and (2) have at least one file, and then returns a set of
    folder paths.
    
    '''
    
    dir_contents = os.listdir(path)
    folder_found = 0
    if len(dir_contents) > 0:
        for temp_object in dir_contents:
            new_path = os.path.join(path, temp_object)
            if os.path.isdir(new_path) == True:
                new_folders = find_terminal_folders(new_path, folders)
                folders = folders.union(new_folders)
                folder_found = 1
        if folder_found == 0:
            folders = folders.union(set([path]))
    else:
        folders = folders.union(set([path]))
        
    return folders
    
def process_new_subjects(batch_size=100,
                         input_bucket_name = 'hbcd-main-study',
                         output_bucket_name = 'hbcd-dicoms-main',
                         input_bucket_config = None,
                         output_bucket_config = None,
                         base_directory = None,
                         symri_container_path = '',
                         layout_path = '',
                         global_path = '',
                         dcm2bids_config_path = '',
                         check_date = True,
                         clean_up = True):
    '''Wrapper function to run symri processing
    
    This function was designed to manage processing of SyMRI QALAS data
    for HBCD. The function assumes that there is some input bucket
    where the only contents ending in *gz are compressed folders
    containing dicoms. It is assumed that the structure is of
    the form /series_number/DICOM/images/dcm. While looking through
    the bucket, the function will only be searching for files that
    (1) have been uploaded more than a day ago, and (2) have not already
    been processed. By default a maximum of 100 gz files will be downloaded,
    but this can be changed by modifying batch_size. The gz files will
    be downloaded to a folder with s3_cache in its name under the base_directory
    specified by the user. Once all gz files are downloaded, the script will go
    one by one through each file and (1) look for dicom folders with 640 images,
    (2) uncompress the files corresponding with those folders, (3) run the most
    recent series through symri's container, and (4) through BIDS conversion.
    After this has been done for all sessions, the data will be uploaded to the
    S3 bucket, and a new tracking log will be created to describe the gz files
    that had been processed. At this point all local files created by processing
    will be deleted leaving only the tracking log and the contents in s3. The next
    time this function is ran, the latest tracking log will be referenced to determine
    what gz files need processing. In the case that this script is stopped during
    processing or the processed data is unable to be uploaded to S3, the tracking
    log will not be updated (except for subjects that were uploaded). Running this
    function again will therefore redo processing for those subjects. If there is
    a gz file that shouldn't be processed, the user could add a new line to the log
    file with the name of the gz file so that this function thinks the gz file has
    already been processed.
    
    Parameters
    ----------
    
    batch_size: int or None
        number of subjects to process, default 100.
        If None, then all subjects will be processed
    input_bucket_name: str
        The name of the existing bucket to grab gz
        files from
    output_bucket_name: str
        The name of the existing bucket to store
        BIDS outputs
    base_directory: str
        The working directory to house all processing
        and tracking logs. A tracking log mush exist at
        /base_directory/tracking_logs/tracking_log*.csv
        to run processing (if none exists, create empty
        csv with columns specified in update_tracking_log
        without any spaces)
    check_date: boolean, default True
        If True, only gz files that have been in the bucket
        for at least one day will be downloaded
    clean_up: boolean, default True
        Remove working directories after results have been
        sent to S3. Only set to False if you want to do
        debugging.
    
    '''

    if os.path.exists(base_directory) == False:
        os.makedirs(base_directory)
    working_dir_base_path = os.path.join(base_directory, 'working_dir')
    if os.path.exists(working_dir_base_path) == False:
        os.makedirs(working_dir_base_path)
    tracking_log_dir = os.path.join(base_directory, 'tracking_logs')
    if os.path.exists(tracking_log_dir) == False:
        os.makedirs(tracking_log_dir)
    
    ############################################################################################################
    ############################################################################################################
    #1. Find which files have already been processed
    #Grab the most recent tracking log
    tracking_logs = glob.glob(os.path.join(base_directory,'tracking_logs/tracking_log*.csv'))
    if len(tracking_logs) == 0:
        template_tracking_log =  os.path.join(Path(inspect.getfile(main)).absolute().parent.resolve(), 'empty_tracking_log.csv')
        shutil.copyfile(template_tracking_log, os.path.join(tracking_log_dir, 'empty_tracking_log.csv'))
        most_recent_tracking_log = os.path.join(tracking_log_dir, 'empty_tracking_log.csv')
    else:
        tracking_logs.sort()
        most_recent_tracking_log = tracking_logs[-1]

    tracking_log = pd.read_csv(most_recent_tracking_log)
    processed_tgz_files = tracking_log['tar_name'].values

    ############################################################################################################
    ############################################################################################################
    #2. Download new tgz files, cross referencing against the files
    # that have already been processed

    date_time = datetime.now().strftime("_date%Y_%m_%d_time%H_%M_%S")
    s3_cache_dir = os.path.join(base_directory, 's3_cache', 's3_cache{}'.format(date_time))
    targz_paths = grab_from_s3(processed_tgz_files, batch_size = batch_size, bucket_name = input_bucket_name, extension = 'gz',
                        output_folder = s3_cache_dir, check_date = check_date, different_config_path = input_bucket_config)
    if len(targz_paths) == 0:
        sys.stdout.write('No new tar files found for processing!\n\n')
        sys.stdout.flush()
        return

    ############################################################################################################
    ############################################################################################################
    #3. Iterate through the gz files, running processing on each one

    output_info = []
    for temp_path in targz_paths:
        sys.stdout.write('Checking if the following archive should be processed: {}\n'.format(temp_path.split('/')[-1]))
        sys.stdout.flush()
        output_info.append(convert_single_tar(temp_path,
                                              working_dir_base_path,
                                              symri_container_path,
                                              layout_path,
                                              global_path,
                                              dcm2bids_config_path))

    ############################################################################################################
    ############################################################################################################
    #4. Upload the converted files to S3

    upload_status = []
    for i in range(len(output_info)):

        #Don't upload anything if the QALAS -> Synthetic Map Conversion failed
        if output_info[i]['symri_conversion_error'] == 1:
            upload_status.append(True)
        
        #Upload if conversion was succesful
        elif output_info[i]['converted_series_name'] != '-1':
            temp = push_to_s3(output_info[i]['bids_path'], output_info[i]['subject_label'], bucket_name = output_bucket_name,
                                prefix = os.path.join('derivatives', output_info[i]['session_label'], 'symri'),
                                different_config_path=output_bucket_config)
            upload_status.append(temp)
            
        #Also don't upload anything if there weren't any QALAS scans found
        else:
            upload_status.append(True)


    #############################################################################################################
    #############################################################################################################
    #5. Create a new log file

    sys.stdout.write('Output Info: {}\n\n'.format(output_info))
    sys.stdout.flush()
    update_tracking_log(most_recent_tracking_log, output_info, upload_status, base_path = os.path.join(base_directory, 'tracking_logs'))
    sys.stdout.write('New tracking log at: ' + most_recent_tracking_log + '\n\n')
    sys.stdout.flush()

    #############################################################################################################
    #############################################################################################################
    #6. Clean up working directories

    if clean_up:
        for i in range(len(output_info)):

            #This will already have been deleted if processing wasnt ran
            if os.path.exists(output_info[i]['working_dir_path']):
                shutil.rmtree(output_info[i]['working_dir_path'])

        shutil.rmtree(s3_cache_dir)
    
    return


def main():
    parser = argparse.ArgumentParser(description="Workflow to Make Relaxometry Maps from QALAS Images Using SyMRI Container.")
    parser.add_argument("ucsd_config_path", help="Path to ucsd s3 config file.")
    parser.add_argument("loris_config_path", help="Path to loris s3 config file.")
    parser.add_argument("symri_container_path", help="Path to SyMRI container.")
    parser.add_argument("base_directory_for_proc", help="Path where processing will occur and where logs will be stored.")
    parser.add_argument('--custom_symri_layout', help="Path to a non-default SyMRI layout file.", type=str)
    parser.add_argument('--custom_global_path', help="Path to a non-default SyMRI global file.", type=str)
    parser.add_argument('--custom_dcm2bids_config_path', help="Path to a non-default dcm2bids config file.", type=str)
    parser.add_argument('--custom_processing_batch_size', help="The number of dicom archives that you want to attempt to process.", type=int, default=20)
    parser.add_argument('--custom_dicom_bucket_name', help="The name of the bucket to grab dicoms from.", type=str, default='midb-hbcd-ucsd-main-pr-dicoms')
    parser.add_argument('--custom_loris_bucket_name', help="The name of the bucket where results will be stored.", type=str, default='midb-hbcd-main-pr')
    parser.add_argument('--clean_up', help="If used, local copies of dicom/niftis generated during processing will be removed.", action='store_true')
    parser.add_argument('--check_date', help="If used, data will only be processed if it is at least ~1 day old.", action='store_true')
    args = parser.parse_args()


    #Set layout file
    if args.custom_symri_layout:
        symri_layout = args.custom_symri_layout
    else:
        symri_layout =  os.path.join(Path(inspect.getfile(main)).absolute().parent.resolve(), 'parametricMaps_Sag_layout.ini')

    #Set global file
    if args.custom_global_path:
        symri_global = args.custom_global_path
    else:
        symri_global =  os.path.join(Path(inspect.getfile(main)).absolute().parent.resolve(), 'global.ini')

    #Set dcm2bids config file
    if args.custom_dcm2bids_config_path:
        dcm2bids_config = args.custom_dcm2bids_config_path
    else:
        dcm2bids_config =  os.path.join(Path(inspect.getfile(main)).absolute().parent.resolve(), 'dcm2bids_config.json')
    
    process_new_subjects(batch_size=args.custom_processing_batch_size,
                         input_bucket_name = args.custom_dicom_bucket_name,
                         output_bucket_name = args.custom_loris_bucket_name,
                         input_bucket_config = args.loris_config_path,
                         output_bucket_config = args.ucsd_config_path,
                         base_directory = args.base_directory_for_proc,
                         symri_container_path = args.symri_container_path,
                         layout_path = symri_layout,
                         global_path = symri_global,
                         dcm2bids_config_path = dcm2bids_config,
                         check_date = args.check_date,
                         clean_up = args.clean_up)


if __name__ == "__main__":
    main()


