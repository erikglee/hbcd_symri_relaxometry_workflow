import boto3
import numpy as np
import os, glob, json, sys, tarfile, shutil
from pydicom import dcmread
from datetime import datetime
from datetime import timezone
from pathlib import Path
import inspect
import argparse


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


def create_page_iterator(bucket = 'hbcd-pilot', prefix = 'derivatives', bids_bucket_config = False):
    '''Utility to create a page iterator for s3 bucket'''
    
    #Grab config path
    if bids_bucket_config == False:
        config_path = ''
    else:
        if type(bids_bucket_config) != str:
            raise NameError('Error: different config path should eithe be string or boolean')
        else:
            config_path = bids_bucket_config
            
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

    del access_key, secret_key, host_base
        
    # Create a reusable Paginator
    paginator = client.get_paginator('list_objects')
        
    # Create a PageIterator from the Paginator
    page_iterator = paginator.paginate(Bucket=bucket, Prefix=prefix)
    
    return page_iterator, client

def get_file_names_with_ending(dicom_bucket_name, dicom_bucket_config, ending = '.json', prefix = ''):
    """Retrieves the names of all JSON files from an S3 bucket.

    Args:
        bucket_name (str): The name of the S3 bucket.

    Returns:
        list: A list of the names of the JSON files in the bucket.
    """

    page_iterator, _ = create_page_iterator(bucket = dicom_bucket_name, prefix = prefix, bids_bucket_config = dicom_bucket_config)

    json_file_names = []
    for page in page_iterator:
        if 'Contents' in page:
            for obj in page['Contents']:
                if obj['Key'].endswith(ending):
                    json_file_names.append(obj['Key'])

    return json_file_names

def download_s3_files_by_name(bids_bucket_config, output_folder, files_to_download, bucket = 'hbcd-pilot'):
    '''Utility to find BIDS subjects in S3 bucket
    
    Parameters
    ----------
    
    bids_bucket_config : str
        This will be used as a config file to identify
        the s3 credentials
    bucket : str, default 'hbcd-pilot'
        The name of the bucket to query for subjects
    prefix : str, default 'assembly_bids'
        The prefix to restrict the files returned by
        the search query (i.e. if data is at
        s3://hbcd-pilot/assembly_bids/sub-1, then
        prefix = 'assembly_bids')
        
    Returns
    -------
    
    s3_subjects : list
        List of "subjects", meaning any instance where
        the path of a folder in the s3 bucket had 'sub-'
        afte the first '/' following the prefix.
        This doesn't mean that the file actually
        corresponds to a BIDS subject, just that it satisfies
        this simple naming pattern.
        
    '''

    # Create a PageIterator    
    _, client = create_page_iterator(bucket = bucket, prefix = '', bids_bucket_config = bids_bucket_config)

    #Iterate through bucket to find potential subjects
    downloaded_files = []
    for temp_file in files_to_download:

        parent_dir = os.path.dirname(os.path.join(output_folder, temp_file))
        if os.path.exists(parent_dir) == False: os.makedirs(parent_dir)
        #print('Bucket Name: {}'.format(bucket))
        #print('File to download: {}'.format(temp_file))
        #print('Output path: {}'.format(os.path.join(output_folder, temp_file)))
        #print(os.path.join(output_folder, temp_file))
        #client.download_file(bucket, temp_dict['Key'], os.path.join(output_folder, temp_dict['Key']))
        #downloaded_files.append(temp_dict['Key'])
        client.download_file(bucket, temp_file, os.path.join(output_folder, temp_file))
        downloaded_files.append(temp_file)
            
    #return list of s3 subjects
    return downloaded_files

def update_best_qalas_info_dict(scan_dict, json_name):
    
    best_qalas_info = {}
    best_qalas_info['QU_motion'] = scan_dict['QU_motion']
    best_qalas_info['aqc_motion'] = scan_dict['aqc_motion']
    best_qalas_info['SeriesInstanceUID'] = scan_dict['SeriesInstanceUID']
    best_qalas_info['StudyInstanceUID'] = scan_dict['StudyInstanceUID']
    best_qalas_info['SubjID'] = scan_dict['SubjID']
    best_qalas_info['json_name'] = json_name
    best_qalas_info['archive_to_download'] = json_name.replace('_mripcqc_info.json', '.tar.gz').split('/')[-1]
    
    return best_qalas_info

def unpack_qalas_from_targz(tar_path, output_path, SeriesInstanceUID = None, StudyInstanceUID = None):
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
    
    #try:    
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

                    #If the user specified a series instance uid, only
                    #use if this matches what is observed in the dicom
                    print('A')
                    if type(SeriesInstanceUID) != type(None):
                        if str(tmp_dcm[0x0020, 0x000e]._value) != SeriesInstanceUID:
                            continue

                    #If the user specified a study instance uid,
                    #only use if this matches what is observed in the dicom
                    if type(StudyInstanceUID) != type(None):
                        if str(tmp_dcm[0x0020, 0x000d]._value) != StudyInstanceUID:
                            continue
                    qalas_folders.append(temp_dir)
                    temp_dict = {}
                    temp_dict['SeriesInstanceUID'] = str(tmp_dcm[0x0020, 0x000e]._value)
                    temp_dict['StudyInstanceUID'] = str(tmp_dcm[0x0020, 0x000d]._value)
                    temp_dict['PatientName'] = str(tmp_dcm[0x0010, 0x0010]._value)
                    temp_dict['PatientID'] = str(tmp_dcm[0x0010, 0x0020]._value)
                    supplemental_infos.append(temp_dict.copy())
                    break
    #except:
    #print('Error: Unable to successfully parse through the dicoms of the current archive. Skipping processing for this archive.')
    #sys.stdout.write('Error: Unable to successfully parse through the dicoms of the current archive. Skipping processing for: {}\n\n'.format(tar_path))
    #sys.stdout.flush()
    #return None, None
        
    return qalas_folders, supplemental_infos

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


def convert_single_tar(qalas_folders, supplemental_infos, qalas_info_dict,
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

    output_info = {'tar_name' : qalas_info_dict['archive_to_download']}
    
    #Create working directory with time stamp to store BIDS output
    tmp_bids_dir = os.path.join(working_dir_base_path, 'tmp_bids')
    if os.path.exists(tmp_bids_dir) == False:
        sys.stdout.write('Making working directory at: {}\n'.format(tmp_bids_dir))
        sys.stdout.flush()
        os.makedirs(tmp_bids_dir)

    #Parse the targz file name to find the subject and session label
    [subject_label, session_label] = targz_to_sub_ses_labels(qalas_info_dict['archive_to_download'])
    
    output_info['bids_path'] = tmp_bids_dir
    date_time = datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S")
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
    if type(qalas_folders) == type(None):
        output_info['converted_series_name'] = str(-1)
        output_info['num_qalas_scans'] = -1
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 1
        shutil.rmtree(output_info['working_dir_path'])
        return output_info
        
    if len(qalas_folders) == 1:
        output_info['converted_series_name'] = qalas_folders
    else:
        raise NameError('Error: Script is trying to unpack data for the following archive but multiple qalas archives were erroneously identified: {}'.format(qalas_info_dict['archive_to_download']))
    
    
    initial_log_path = os.path.join(working_dir_base_path, 'symri_container.log')
    
    #Specify the base commands for qalas and dcm2bids that will be configured for each qalas folder
    qalas_base_command = 'singularity run -B {global_path}:/opt/symri/bin/global.ini -B {layout_path}:/layout_path -B {qalas_folder}:/input -B {dcm_maps_path}:/output {container_path} --batch-mode --dataset /input --output /output --layout /layout_path --force-anatomy "infant brain" &> {log_path}'
    dcm2bids_base_command = '{dcm2bids_executable_path} -d {dcm_maps_path} --participant {subject_label} --session {session_label} -c {dcm2bids_config_path} -o {tmp_bids_dir}'
    dcm2bids_executable_path = 'dcm2bids'

    #Run symri container
    dcm_maps_path = os.path.join(working_dir_base_path, qalas_folders[0] + '_qalas_derived_dcm_maps')
    os.makedirs(dcm_maps_path)
    symri_container_command = qalas_base_command.format(global_path=global_path, layout_path=layout_path, container_path=container_path, qalas_folder=qalas_folders[0], dcm_maps_path=dcm_maps_path, log_path=initial_log_path)
    os.system(symri_container_command)
    sys.stdout.write('SyMRI Command:\n')
    sys.stdout.write(symri_container_command + '\n\n')
        
    #Run dcm2bids conversion
    output_info['num_niftis_generated'] = 0
    print(dcm_maps_path)
    sys.stdout.write('dcm_maps_path: {}\n'.format(dcm_maps_path))
    sys.stdout.flush()
    if len(glob.glob(os.path.join(dcm_maps_path, '*'))) == 6:

        dcm2bids_command = dcm2bids_base_command.format(dcm2bids_executable_path=dcm2bids_executable_path, dcm_maps_path=dcm_maps_path, dcm2bids_config_path=dcm2bids_config_path, tmp_bids_dir=tmp_bids_dir, subject_label=subject_label, session_label=session_label)
        os.system(dcm2bids_command)
        shutil.copyfile(initial_log_path, os.path.join(tmp_bids_dir, 'sub-' + subject_label, 'ses-' + session_label, 'anat', 'sub-{}_ses-{}_acq-QALAS_desc-SymriContainer.log'.format(output_info['subject_label'], output_info['session_label'])))

        output_info['symri_conversion_error'] = 0
        output_info['num_niftis_generated'] += len(glob.glob(os.path.join(tmp_bids_dir, 'sub*','ses*','anat','*.nii.gz')))

        #Update the metadata in the json files
        json_files = glob.glob(os.path.join(tmp_bids_dir, 'sub-' + subject_label, 'ses-' + session_label, 'anat', '*json'))
        for temp_json_file in json_files:
            update_bids_json(temp_json_file, supplemental_info=supplemental_infos[0])

        sys.stdout.write('dcm2bids Command:\n')
        sys.stdout.write(dcm2bids_command + '\n\n')
        sys.stdout.flush()



    else:
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 1
        sys.stdout.write('SyMRI Container was unable to convert QALAS dicoms to Synthetic Maps\n')
        sys.stdout.flush()
    
    return output_info


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

    #sys.stdout.write(base_bids_dir + '\n')
    #sys.stdout.write(subject_label + '\n')
    #sys.stdout.write(bucket_name + '\n')
    #sys.stdout.write(prefix + '\n')
    #sys.stdout.write(different_config_path + '\n')
    #sys.stdout.flush()
    

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

    #sys.stdout.write('Client Info: {}\n'.format(client))
    #sys.stdout.flush()
    
    try:
        sys.stdout.write('Uploading data for sub-{}\n'.format(subject_label))
        sys.stdout.flush()
        os.chdir(base_bids_dir)
        files = glob.glob('sub*/ses*/anat/*')
        sys.stdout.write('Files to upload: {}\n'.format(files))
        for temp_file in files:
            response = client.upload_file(temp_file, bucket_name, os.path.join(prefix, temp_file))
    except:
        sys.stdout.write('Either part or all of uploading Failed for sub-{}, the corresponding session will not be save to the log file so that processing will run again later.\n'.format(subject_label))
        sys.stdout.flush()
        return False

    return True
    

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
    parser.add_argument('--keep_work_dirs', help="If used, local copies of dicom/niftis generated will be saved. Deletion may be required for subsequent processing in specific cases.", action='store_true')
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

    #First grab all JSONS and group them by session
    files = get_file_names_with_ending(args.custom_dicom_bucket_name, args.ucsd_config_path, ending = '_mripcqc_info.json', prefix = '')
    jsons_dict = {}
    for temp_file in files:
        file_split = '_'.join(temp_file.split('_')[0:3])
        if file_split not in jsons_dict.keys():
            jsons_dict[file_split] = [temp_file]
        else:
            jsons_dict[file_split].append(temp_file)


    #See if a tracking log exists and if so load it. Tracking
    #log should be in the form of a dictionary of dictionaries.
    logs_dir = os.path.join(args.base_directory_for_proc, 'logs')
    if os.path.exists(logs_dir) == False:
        os.makedirs(logs_dir)
    log_files = glob.glob(os.path.join(logs_dir, 'tracking_log*.json'))
    log_files.sort()
    if len(log_files) == 0:
        tracking_log_data = {}
    else:
        with open(log_files[-1], 'r') as f:
            tracking_log_data = json.load(f)
            
    #Second, figure out which jsons (1) represent sessions that need to be
    #processed, (2) sessions that may need to be modified, and (3) sessions
    #that are already set.
    sessions_to_evaluate = {}
    potentially_need_reprocessing = []
    batch_limit = args.custom_processing_batch_size
    num_to_process = 0
    for temp_session in jsons_dict.keys():
        if temp_session in tracking_log_data.keys():
            already_processed_archives = tracking_log_data[temp_session]['jsons_tested']
            num_same = 0
            for temp_json in jsons_dict[temp_session]:
                tar_name = temp_json.replace('_mripcqc_info.json', '.tar.gz')
                if tar_name in already_processed_archives:
                    num_same += 1
            if num_same == len(already_processed_archives):
                #No need to do anything to update subject
                continue
            else:
                potentially_need_reprocessing.append(temp_session)
        else:
            if num_to_process < batch_limit:
                
                sessions_to_evaluate[temp_session] = jsons_dict[temp_session]
                num_to_process += 1
                
        if num_to_process == batch_limit:
            break
                
    #Third iterate through all sessions that need to be processed. #First,
    #download the QC JSONS for the sessions. 
    archives_with_qalas_to_process = {}
    working_base_dir = os.path.join(args.base_directory_for_proc, 'working_dir')
    for temp_session in sessions_to_evaluate.keys():
        
        session_base_dir = os.path.join(working_base_dir, temp_session)
        dicom_base_dir = os.path.join(session_base_dir, 'dicoms')
        jsons_base_dir = os.path.join(session_base_dir, 'qc_jsons') 
        if os.path.exists(jsons_base_dir) == False:
            os.makedirs(jsons_base_dir)
        download_s3_files_by_name(args.ucsd_config_path, jsons_base_dir, jsons_dict[temp_session], bucket = args.custom_dicom_bucket_name)
        downloaded_jsons = glob.glob(os.path.join(jsons_base_dir, '*.json'))
        
        best_qalas_info = None
        best_qalas_qu_score = np.inf
        best_qalas_aqc_score = np.inf
        for temp_json in downloaded_jsons:
            with open(temp_json, 'r') as f:
                content = json.load(f)[0]
            if type(content) == dict:
                content = [content]
            for temp_scan in content:
                if temp_scan['SeriesType'] == 'qMRI':
                    if (type(temp_scan['QU_motion']) == type(None)) or (type(temp_scan['aqc_motion']) == type(None)) or (type(temp_scan['Completed']) == type(None)):
                        print('Either QU_motion, aqc_motion, or Completed status is missing for one scan within {}'.format(temp_json))
                        continue
                    if temp_scan['Completed'] == 1:
                        if temp_scan['QU_motion'] < best_qalas_qu_score:
                            best_qalas_qu_score = temp_scan['QU_motion']
                            best_qalas_aqc_score = temp_scan['aqc_motion']
                            best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)
                        elif temp_scan['QU_motion'] == best_qalas_qu_score:
                            if temp_scan['aqc_motion'] < best_qalas_aqc_score:
                                best_qalas_qu_score = temp_scan['QU_motion']
                                best_qalas_aqc_score = temp_scan['aqc_motion']
                                best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)
                                
        if type(best_qalas_info) != type(None):
            
            #Download the tar archive that has the best qalas scan
            best_qalas_info['jsons_tested'] = sessions_to_evaluate[temp_session]
            if os.path.exists(dicom_base_dir) == False:
                os.makedirs(dicom_base_dir)
            file_names = download_s3_files_by_name(args.ucsd_config_path, dicom_base_dir, [best_qalas_info['archive_to_download']], bucket = args.custom_dicom_bucket_name)
            archives_with_qalas_to_process[temp_session] = best_qalas_info
            
            #Unpack the archive and identify which folder has the qalas scan
            unpack_dir = os.path.join(session_base_dir, 'dicoms', 'unpacked_archive')
            if os.path.exists(unpack_dir) == False:
                os.makedirs(unpack_dir)
            tar_path = os.path.join(dicom_base_dir, archives_with_qalas_to_process[temp_session]['archive_to_download'])
            qalas_folders, supplemental_infos = unpack_qalas_from_targz(tar_path, unpack_dir, SeriesInstanceUID = archives_with_qalas_to_process[temp_session]['SeriesInstanceUID'], StudyInstanceUID = archives_with_qalas_to_process[temp_session]['StudyInstanceUID'])    
            
            print(qalas_folders)
            print(supplemental_infos)
            
            #Run SyMRI to create relaxometry maps and convert to BIDS
            output_info = convert_single_tar(qalas_folders, supplemental_infos, best_qalas_info,
                                            session_base_dir,
                                            args.symri_container_path,
                                            symri_layout,
                                            symri_global,
                                            dcm2bids_config)
            
            #Send the results to s3
            status = push_to_s3(output_info['bids_path'], output_info['subject_label'], bucket_name = args.custom_loris_bucket_name,
                                    prefix = os.path.join('derivatives', 'ses-' + output_info['session_label'], 'symri'),
                                    different_config_path=args.loris_config_path)     
            if status == False:
                raise ValueError('Error: Pushing data to S3 was unsuccessful.')
            
            #After BIDS conversion + uploading to S3, update the
            #tracking dictionary with data for this subject.
            tracking_log_data[temp_session] = best_qalas_info
            
        
        else:
            
            #Add the jsons for this session to the json log so they
            #can be ignored in further processing.
            best_qalas_info = {}
            best_qalas_info['jsons_tested'] = jsons_dict[temp_session]
            tracking_log_data[temp_session] = best_qalas_info

        #Clean up the session working directory
        if args.keep_work_dirs == False:
            shutil.rmtree(session_base_dir)
            
    #Save the tracking log
    date_time = datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S")
    new_path = os.path.join(logs_dir, 'tracking_log_{}.json'.format(date_time))
    with open(new_path, 'w') as f:
        f.write(json.dumps(tracking_log_data, indent = 5))

    return


if __name__ == "__main__":
    main()

