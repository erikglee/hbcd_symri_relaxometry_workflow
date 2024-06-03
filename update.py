import boto3
import numpy as np
import os, glob, json, sys, tarfile, shutil
from pydicom import dcmread
from datetime import datetime
from datetime import timezone
from pathlib import Path
import inspect
import argparse
import time


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
    details = {"ImageDescription" : "This is a synthetic image derived from a QALAS scan distributed by SyMRI.\\nQunatitative T1, T2, and PD values are estimated from the QALAS scan using numerical algorithms provided by SyMRI. B1 values are estimated by SyMRI.",
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
        s3://bucket_name/assembly_bids/sub-1, then
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
    try:
        downloaded_files = []
        for temp_file in files_to_download:

            parent_dir = os.path.dirname(os.path.join(output_folder, temp_file))
            if os.path.exists(parent_dir) == False: os.makedirs(parent_dir)
            client.download_file(bucket, temp_file, os.path.join(output_folder, temp_file))
            downloaded_files.append(temp_file)
    except:
        sys.stdout.write('Error: Unable to download the following file: {}\n'.format(temp_file))
        sys.stdout.flush()
        return None
            
    #return list of s3 subjects
    return downloaded_files

def update_best_qalas_info_dict(scan_dict, json_name):
    
    best_qalas_info = {}
    best_qalas_info['QU_motion'] = scan_dict['QU_motion']
    best_qalas_info['aqc_motion'] = scan_dict['aqc_motion']
    best_qalas_info['SeriesInstanceUID'] = scan_dict['SeriesInstanceUID']
    best_qalas_info['StudyInstanceUID'] = scan_dict['StudyInstanceUID']
    best_qalas_info['SubjID'] = scan_dict['SubjID']
    best_qalas_info['json_for_unpacked_archive'] = json_name
    best_qalas_info['archive_to_download'] = json_name.replace('_mripcqc_info.json', '.tar.gz').split('/')[-1]
    
    return best_qalas_info

def qalas_selection_with_qu_motion(downloaded_jsons):
    best_qalas_info = None
    best_qalas_qu_score = np.inf
    best_qalas_aqc_score = np.inf
    best_qalas_snr = -1*np.inf
    for temp_json in downloaded_jsons:
        print('         For loop on jsons')
        with open(temp_json, 'r') as f:
            content = json.load(f)[0]
        if type(content) == dict:
            content = [content]
        for temp_scan in content:
            if temp_scan['SeriesType'] == 'qMRI':
                sys.stdout.write('         qMRI present in current json: {}\n'.format(temp_json))
                sys.stdout.write('            QU: {}\n'.format(temp_scan['QU_motion']))
                sys.stdout.write('            Compliant: {}\n'.format(temp_scan['HBCD_compliant']))
                sys.stdout.write('            AQ: {}\n'.format(temp_scan['aqc_motion']))
                sys.stdout.write('            brain_SNR: {}\n'.format(temp_scan['brain_SNR']))
                sys.stdout.flush()
                #Be sure that at least the QU_motion, aqc_motion and HBCD_compliant fields are there
                if (type(temp_scan['QU_motion']) == type(None)) or (type(temp_scan['aqc_motion']) == type(None)) or (type(temp_scan['HBCD_compliant']) == type(None)) or (type(temp_scan['brain_SNR']) == type(None)):
                    sys.stdout.write('   Processing will either be attempted later or without QU_motion: Either QU_motion, aqc_motion, or HBCD_compliant status is missing for one scan within {}'.format(temp_json))
                    sys.stdout.flush()
                    return None, True
                try:
                    if temp_scan['HBCD_compliant'] == 'Yes':
                        if temp_scan['QU_motion'] < best_qalas_qu_score:
                            best_qalas_qu_score = temp_scan['QU_motion']
                            best_qalas_aqc_score = temp_scan['aqc_motion']
                            best_qalas_snr = temp_scan['brain_SNR']
                            best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)
                        elif temp_scan['QU_motion'] == best_qalas_qu_score:
                            if temp_scan['aqc_motion'] < best_qalas_aqc_score:
                                best_qalas_qu_score = temp_scan['QU_motion']
                                best_qalas_aqc_score = temp_scan['aqc_motion']
                                best_qalas_snr = temp_scan['brain_SNR']
                                best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)
                            elif temp_scan['aqc_motion'] == best_qalas_aqc_score:
                                if temp_scan['brain_SNR'] > best_qalas_snr:
                                    best_qalas_qu_score = temp_scan['QU_motion']
                                    best_qalas_aqc_score = temp_scan['aqc_motion']
                                    best_qalas_snr = temp_scan['brain_SNR']
                                    best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)    
                except Exception as err:
                    print(f"Found the following: {err=}, {type(err)=}")
                    print('Identified while trying to evaluate the following JSON. Likely one or more QC value is missing {}'.format(temp_json.split('/')[-1]))
                    return None, True
                
                            
    return best_qalas_info, False


def qalas_selection_without_qu_motion(downloaded_jsons):
    best_qalas_info = None
    best_qalas_aqc_score = np.inf
    best_qalas_snr = -1*np.inf
    for temp_json in downloaded_jsons:
        with open(temp_json, 'r') as f:
            content = json.load(f)[0]
        if type(content) == dict:
            content = [content]
        for temp_scan in content:
            if temp_scan['SeriesType'] == 'qMRI':
                #Be sure that at least the QU_motion, aqc_motion and HBCD_compliant fields are there
                if (type(temp_scan['aqc_motion']) == type(None)) or (type(temp_scan['HBCD_compliant']) == type(None)) or (type(temp_scan['brain_SNR']) == type(None)):
                    sys.stdout.write('   Processing will either be attempted later or without QU_motion: Either QU_motion, aqc_motion, or HBCD_compliant status is missing for one scan within {}'.format(temp_json))
                    sys.stdout.flush()
                    return None, True
                try:
                    sys.stdout.flush()
                    if temp_scan['HBCD_compliant'] == 'Yes':
                        if temp_scan['aqc_motion'] < best_qalas_aqc_score:
                            best_qalas_aqc_score = temp_scan['aqc_motion']
                            best_qalas_snr = temp_scan['brain_SNR']
                            best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)
                        elif temp_scan['aqc_motion'] == best_qalas_aqc_score:
                            if temp_scan['brain_SNR'] > best_qalas_snr:
                                best_qalas_aqc_score = temp_scan['aqc_motion']
                                best_qalas_snr = temp_scan['brain_SNR']
                                best_qalas_info = update_best_qalas_info_dict(temp_scan, temp_json)    
                except Exception as err:
                    print(f"Found the following: {err=}, {type(err)=}")
                    print('Identified while trying to evaluate the following JSON. Likely one or more QC value is missing {}'.format(temp_json.split('/')[-1]))
                    return None, True
                
                            
    return best_qalas_info, False

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

    If the user specifies a SeriesInstanceUID, or 
    StudyInstanceUID, then only DICOMS that match
    these values will be returned.
    
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
                    temp_dict['PatientAge'] = str(tmp_dcm[0x0010, 0x1010]._value)
                    temp_dict['PatientSex'] = str(tmp_dcm[0x0010, 0x0040]._value)
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
        if (len(split_name[0]) == 9) and (len(split_name[1]) == 6) and (len(split_name[2]) == 3) and (split_name[2][0] == 'V'):
            if len(split_name) == 3:
                sub_label = split_name[1]
                ses_label = split_name[2].split('.')[0]
            else:
                sub_label = split_name[1]
                ses_label = split_name[2]
        else:
            sys.stdout.write('   Unable to parse subject/session labels from the following archive name: {}\n'.format(targz_name.split('/')[-1]))
            sys.stdout.write('   Scripts were expecting 9 char PSCID, 6 char DCCID, and 3 char session label starting with V for visit, all seperated by underscores (with optional additional text to follow).\n')
            sys.stdout.flush()
            
    return sub_label, ses_label


def convert_single_tar(qalas_folders, supplemental_infos, qalas_info_dict,
                       working_dir_base_path,
                       container_path,
                       layout_path,
                       global_path,
                       dcm2bids_config_path):
    
    '''Function to unpack qalas from dcm, convert to maps then bids
    
    '''

    output_info = {'tar_name' : qalas_info_dict['archive_to_download']}
    
    #Create working directory with time stamp to store BIDS output
    tmp_bids_dir = os.path.join(working_dir_base_path, 'tmp_bids')
    if os.path.exists(tmp_bids_dir) == False:
        sys.stdout.write('   Making working directory at: {}\n'.format(tmp_bids_dir))
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
        return output_info

    #Unpack ONLY the qalas folders from targz
    if type(qalas_folders) == type(None):
        output_info['converted_series_name'] = str(-1)
        output_info['num_qalas_scans'] = -1
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 1
        return output_info
        
    if len(qalas_folders) == 1:
        output_info['converted_series_name'] = qalas_folders
    else:
        raise NameError('Error: Script is trying to unpack data for the following archive but multiple qalas archives were erroneously identified: {}'.format(qalas_info_dict['archive_to_download']))
    
    #Update the supplmental info dict with QC values.
    supplemental_infos[0]['QU_motion'] = qalas_info_dict['QU_motion']
    supplemental_infos[0]['aqc_motion'] = qalas_info_dict['aqc_motion']
    
    initial_log_path = os.path.join(working_dir_base_path, 'symri_container.log')
    
    #Specify the base commands for qalas and dcm2bids that will be configured for each qalas folder
    qalas_base_command = 'singularity run -B {global_path}:/opt/symri/bin/global.ini -B {layout_path}:/layout_path -B {qalas_folder}:/input -B {dcm_maps_path}:/output {container_path} --batch-mode --dataset /input --output /output --layout /layout_path --force-anatomy "infant brain" &> {log_path}'
    dcm2bids_base_command = '{dcm2bids_executable_path} -d {dcm_maps_path} --participant {subject_label} --session {session_label} -c {dcm2bids_config_path} -o {tmp_bids_dir}'
    dcm2bids_executable_path = 'dcm2bids'

    #Run symri container
    dcm_maps_path = os.path.join(working_dir_base_path, 'dcm_maps_dir')
    os.makedirs(dcm_maps_path)
    symri_container_command = qalas_base_command.format(global_path=global_path, layout_path=layout_path, container_path=container_path, qalas_folder=qalas_folders[0], dcm_maps_path=dcm_maps_path, log_path=initial_log_path)
    os.system(symri_container_command)
    sys.stdout.write('   SyMRI Command:\n')
    sys.stdout.write('      ' + symri_container_command + '\n')
        
    #Run dcm2bids conversion
    output_info['num_niftis_generated'] = 0
    output_info['dcm2bids_conversion_error'] = 0
    sys.stdout.write('   dcm_maps_path: {}\n'.format(dcm_maps_path))
    sys.stdout.flush()
    time.sleep(10)
    if len(glob.glob(os.path.join(dcm_maps_path, '*'))) == 7:

        try:
            dcm2bids_command = dcm2bids_base_command.format(dcm2bids_executable_path=dcm2bids_executable_path, dcm_maps_path=dcm_maps_path, dcm2bids_config_path=dcm2bids_config_path, tmp_bids_dir=tmp_bids_dir, subject_label=subject_label, session_label=session_label)
            os.system(dcm2bids_command)
            sys.stdout.write('   dcm2bids Command:\n')
            sys.stdout.write('      ' + dcm2bids_command + '\n')
            shutil.copyfile(initial_log_path, os.path.join(tmp_bids_dir, 'sub-' + subject_label, 'ses-' + session_label, 'anat', 'sub-{}_ses-{}_acq-QALAS_desc-SymriContainer.log'.format(output_info['subject_label'], output_info['session_label'])))
        except:
            sys.stdout.write('   dcm2bids conversion failed for the following archive: {}\n'.format(qalas_info_dict['archive_to_download']))
            sys.stdout.flush()
            output_info['num_niftis_generated'] += len(glob.glob(os.path.join(tmp_bids_dir, 'sub*','ses*','anat','*.nii.gz')))
            output_info['symri_conversion_error'] = 0
            output_info['dcm2bids_conversion_error'] = 1
            return output_info

        output_info['symri_conversion_error'] = 0
        output_info['num_niftis_generated'] += len(glob.glob(os.path.join(tmp_bids_dir, 'sub*','ses*','anat','*.nii.gz')))

        #Update the metadata in the json files
        json_files = glob.glob(os.path.join(tmp_bids_dir, 'sub-' + subject_label, 'ses-' + session_label, 'anat', '*json'))
        for temp_json_file in json_files:
            update_bids_json(temp_json_file, supplemental_info=supplemental_infos[0])

        sys.stdout.write('   dcm2bids Command:\n')
        sys.stdout.write('      ' + dcm2bids_command + '\n\n')
        sys.stdout.flush()



    else:
        output_info['num_niftis_generated'] = 0
        output_info['symri_conversion_error'] = 1
        sys.stdout.write('   SyMRI Container was unable to convert QALAS dicoms to Synthetic Maps\n')
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
        sys.stdout.write('   Files to upload: {}\n'.format(files))
        for temp_file in files:
            response = client.upload_file(temp_file, bucket_name, os.path.join(prefix, temp_file))
    except:
        sys.stdout.write('   Either part or all of uploading Failed for sub-{}, the corresponding session will not be save to the log file so that processing will run again later.\n'.format(subject_label))
        sys.stdout.flush()
        return False

    return True
    

def main():

    '''Converts QALAS DICOM archives to BIDS derivatives relaxometry maps.

    This script will download QALAS DICOM archives from an S3 bucket, unpack
    them, convert them to relaxometry maps using the SyMRI container, and then
    convert the resulting DICOMs to BIDS derivatives. The script will then
    upload the BIDS derivatives to a different S3 bucket.

    Steps of processing are roughly as follows:
    1. Identify all json files in the dicom bucket that have the ending
         '_mripcqc_info.json' and group them by subject/session. The
         json files have information about the quality of the DICOMs within
         a given archive. If a subject/session combo has multiple archives,
         and/or multiple QALAS scans, the best QALAS scan will be kept for
         processing based on QC info provided by the json file.
    2. If the file reproc_log.json exists in the base directory, then
        load it. This file contains a list of subjects/sessions that could need
        to be reprocessed based on new data. If a subject appears in the
        any portion of the json aside from "to_reprocess", any new or existing
        archives for the session will be ignored. If the user has added a 
        subject/session to the to_reprocess portion of the json, then the code will attempt
        to overwrite any existing processing for that subject/session. This
        may result in either the same or different processing results depending
        on the new data that is available for the subject. If the reproc_log.json
        file does not exist, then the code will create it.
    3. Iterate through all subjects/sessions that need to be processed. For subjects/
        sessions that have a QALAS scan and are ready to be processed, download the
        archive containing the best QALAS scan, unpack the archive, and identify
        the dicom folder that contains the QALAS scan. The best QALAS is determined
        based on the following ranking heirarchy (1) HBCD_Compliant, (2) QU_motion,
        (3) aqc_motion, (4) brain_SNR. If QU_motion is missing from one or more QALAS
        scans, then QU_motion will be excluded from the file comparison procedure. If
        one or more value (besides QU_motion) is missing for one or more QALAS scans,
        processing will not occur.
    4. Run the SyMRI container on the QALAS dicom folder to generate relaxometry maps.
    5. Run dcm2bids on the resulting dicoms to convert them to BIDS derivatives.
        Also copy over a copy of the log file generated by SyMRI during processing. Further,
        add some custom fields to the metadata of the json files generated by dcm2bids that
        will be accompanying the BIDS derivatives.
    6. Upload the BIDS derivatives to a different S3 bucket.
    7. Update the tracking log with information about the processing that was just completed.
    8. Update the reproc_log.json file with information about subjects/sessions that may need
        to be reprocessed in the future. If you want to reprocess a subject, then add the subject
        to the to_reprocess portion of the json.

    If processing is interrupted after partial processing has occurred, the derivatives may still be
    uploaded to S3 despite the tracking logs/jsons not being updated. Upon reprocessing, data from these
    runs will be overwritten which should not cause any issues. If the data under
    {base_directory_for_proc}/work_dir is not properly deleted because of an error, this
    may require the user to delete any contents of the working directory. As long as the log
    files remain intact, processing should be able to pick up in the usual fashion following the deletion.
    
    Parameters
    ----------
    ucsd_config_path : string
        Path to s3 config file for s3 account that has bucket with dicom
        archives (i.e. .tar.gz files)
    loris_config_path : string
        Path to s3 config file for s3 account that has bucket where BIDS
        derivatives will be stored.
    symri_container_path : string
        Path to SyMRI container that will be used for calculating relaxometry
        maps. Container should already have license file set up.
    base_directory_for_proc : string
        Path to directory where processing will occur and where logs will be
        stored. If you change this directory on subsequent calls to this function,
        then processing will start over from scratch.
    custom_symri_layout : string, default None
        Path to a non-default SyMRI layout file.
    custom_global_path : string, default None
        Path to a non-default SyMRI global file.
    custom_dcm2bids_config_path : string, default None
        Path to a non-default dcm2bids config file.
    custom_processing_batch_size : int, default 20
        The number of unique subject/session combos that you want to attempt to process.
    custom_dicom_bucket_name : string, default 'midb-hbcd-ucsd-main-pr-dicoms'
        The name of the bucket to grab dicom archives and json files from. Assumes
        that the each json file/dicom archive is formatted using the format from UCSD
        for the HBCD study.
    custom_loris_bucket_name : string, default 'midb-hbcd-main-pr'
        The name of the bucket where results will be stored. Results will get stored under
        the following prefix: 'derivatives/ses-<session_label>/symri' within the bucket.
    keep_work_dirs : bool, default False
        If used, local copies of dicom/niftis generated will be saved. Deletion may be required
        for subsequent processing in specific cases.
    
    '''

    parser = argparse.ArgumentParser(description="Workflow to Make Relaxometry Maps from QALAS Images Using SyMRI Container.")
    parser.add_argument("ucsd_config_path", help="Path to ucsd s3 config file.")
    parser.add_argument("loris_config_path", help="Path to loris s3 config file.")
    parser.add_argument("symri_container_path", help="Path to SyMRI container.")
    parser.add_argument("base_directory_for_proc", help="Path where processing will occur and where logs will be stored.")
    parser.add_argument('--custom_symri_layout', help="Path to a non-default SyMRI layout file.", type=str)
    parser.add_argument('--custom_global_path', help="Path to a non-default SyMRI global file.", type=str)
    parser.add_argument('--custom_dcm2bids_config_path', help="Path to a non-default dcm2bids config file.", type=str)
    parser.add_argument('--custom_processing_batch_size', help="The number of dicom archives that you want to attempt to process.", type=int, default=20)
    parser.add_argument('--custom_dicom_bucket_name', help="The name of the bucket to grab dicoms from.", type=str, default='hbcd-dicoms-main-study')
    parser.add_argument('--custom_loris_bucket_name', help="The name of the bucket where results will be stored.", type=str, default='midb-hbcd-main-pr')
    parser.add_argument('--keep_work_dirs', help="If used, local copies of dicom/niftis generated will be saved. Deletion may be required for subsequent processing in specific cases.", action='store_true')
    parser.add_argument('--dicom_prefix', help="Prefix to include in front of subject path when searching for dicom archives.", type=str, default='')
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
    files = get_file_names_with_ending(args.custom_dicom_bucket_name, args.ucsd_config_path, ending = '_mripcqc_info.json', prefix = args.dicom_prefix)
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


    #Also look around for a log file in the base dir named reproc_log.json
    #that contains both a list of subjects that could need to be reprocessed
    #based on new data and subjects that the user has deemed ready for reprocessing.
    reproc_log_path = os.path.join(args.base_directory_for_proc, 'reproc_log.json')
    if os.path.exists(reproc_log_path) == True:
        with open(reproc_log_path, 'r') as f:
            reproc_log_data = json.load(f)
        to_reprocess = reproc_log_data['to_reprocess']
        new_archive_available = reproc_log_data['new_archive_available']
        qalas_in_qc_but_not_archive = reproc_log_data['qalas_in_qc_but_not_archive']
        no_niftis_to_upload = reproc_log_data['no_niftis_to_upload']
        missing_archive = reproc_log_data['missing_archive']
    else:
        reproc_log_data = {}
        to_reprocess = []
        new_archive_available = []
        qalas_in_qc_but_not_archive = []
        no_niftis_to_upload = []
        missing_archive = []

    reprocess_attempted = []
    sessions_to_skip = new_archive_available + qalas_in_qc_but_not_archive + no_niftis_to_upload + missing_archive
    for temp_session in to_reprocess:
        if temp_session in sessions_to_skip:
            sessions_to_skip.remove(temp_session)
            
    #Second, figure out which jsons (1) represent sessions that need to be
    #processed, (2) sessions that may need to be modified, and (3) sessions
    #that are already set.
    sessions_to_evaluate = {}
    batch_limit = args.custom_processing_batch_size
    num_to_process = 0
    for temp_session in jsons_dict.keys():
        #Dont process the subject if there was previously something wrong
        #with processing and new processing wasnt explicitly specified by
        #the to_reprocess field.
        if temp_session in sessions_to_skip:
            sys.stdout.write('Skipping session due to status in reproc_log.json file: {}\n'.format(temp_session))
            sys.stdout.flush()
            continue

        #If the session has already been processed, then check to see if
        #any new data has come in since then.
        if temp_session in tracking_log_data.keys():
            already_processed_archives = tracking_log_data[temp_session]['jsons_tested']
            num_same = 0
            for temp_json in jsons_dict[temp_session]:
                if temp_json in already_processed_archives:
                    num_same += 1
            if temp_session in to_reprocess:
                sessions_to_evaluate[temp_session] = jsons_dict[temp_session]
                reprocess_attempted.append(temp_session)
                num_to_process += 1
            elif num_same == len(already_processed_archives):
                #No need to do anything to update subject
                continue
            else:
                new_archive_available.append(temp_session)
        
        #If the session has not been processed, then add it to the list
        else:
            if num_to_process < batch_limit:
                
                sessions_to_evaluate[temp_session] = jsons_dict[temp_session]
                num_to_process += 1
                
        if num_to_process == batch_limit:
            break

    sys.stdout.write('Processing scripts found {} session/subject combonations as candidates for new and/or updated processing.\n'.format(num_to_process))
    sys.stdout.flush()
                
    #Third iterate through all sessions that need to be processed. #First,
    #download the QC JSONS for the sessions. 
    archives_with_qalas_to_process = {}
    working_base_dir = os.path.join(args.base_directory_for_proc, 'working_dir')
    for temp_session in sessions_to_evaluate.keys():

        sys.stdout.write('Attempting processing for: {}\n'.format(temp_session))
        sys.stdout.flush()
        
        session_base_dir = os.path.join(working_base_dir, temp_session)
        dicom_base_dir = os.path.join(session_base_dir, 'dicoms')
        jsons_base_dir = os.path.join(session_base_dir, 'qc_jsons') 
        if os.path.exists(jsons_base_dir) == False:
            os.makedirs(jsons_base_dir)
        downloaded_files = download_s3_files_by_name(args.ucsd_config_path, jsons_base_dir, jsons_dict[temp_session], bucket = args.custom_dicom_bucket_name)
        if type(downloaded_files) == type(None):
            raise ValueError('Error: Unable to download jsons for current subject/session: {}. This should probably never happen if scripts are configured correctly.'.format(temp_session))
        downloaded_jsons = glob.glob(os.path.join(jsons_base_dir, '*.json'))

        sys.stdout.write('   Grabbed the following jsons for current subject/session: {}\n'.format(jsons_dict[temp_session]))
        sys.stdout.flush()
        
        #Evaluate the jsons to find the best QALAS scan
        print('   First evaluating runs based on HBCD_Compliant, QU_motion, aqc_motion, and brain_SNR.')
        best_qalas_info, qc_or_other_missing = qalas_selection_with_qu_motion(downloaded_jsons)
        if type(best_qalas_info) == type(None):
            print('   Instead evaluating runs based only on HBCD_Compliant, aqc_motion, and brain_SNR.')
            best_qalas_info, qc_or_other_missing = qalas_selection_without_qu_motion(downloaded_jsons)

                                
        if type(best_qalas_info) != type(None):

            sys.stdout.write('   Identified QALAS scan with the following info that will be targetted for processing: {}\n'.format(best_qalas_info))
            sys.stdout.flush()
            
            #Download the tar archive that has the best qalas scan
            best_qalas_info['jsons_tested'] = sessions_to_evaluate[temp_session]
            if os.path.exists(dicom_base_dir) == False:
                os.makedirs(dicom_base_dir)
            file_names = download_s3_files_by_name(args.ucsd_config_path, dicom_base_dir, [best_qalas_info['archive_to_download']], bucket = args.custom_dicom_bucket_name)
            if type(file_names) == type(None):
                sys.stdout.write('   Unable to download dicoms for current subject/session: {}. Skipping subject archive and adding them to reprocessing log.\n'.format(temp_session))
                sys.stdout.flush()
                missing_archive.append(temp_session)
                continue
            archives_with_qalas_to_process[temp_session] = best_qalas_info
            
            sys.stdout.write('   The following files have been downloaded: {}\n'.format(file_names))
            sys.stdout.flush()

            #Unpack the archive and identify which folder has the qalas scan
            unpack_dir = os.path.join(session_base_dir, 'dicoms', 'unpacked_archive')
            if os.path.exists(unpack_dir) == False:
                os.makedirs(unpack_dir)
            tar_path = os.path.join(dicom_base_dir, archives_with_qalas_to_process[temp_session]['archive_to_download'])
            qalas_folders, supplemental_infos = unpack_qalas_from_targz(tar_path, unpack_dir, SeriesInstanceUID = archives_with_qalas_to_process[temp_session]['SeriesInstanceUID'], StudyInstanceUID = archives_with_qalas_to_process[temp_session]['StudyInstanceUID'])    
            if len(qalas_folders) == 0:
                sys.stdout.write('   Found QALAS in QC file but unable to find QALAS scan within archive. Skipping subject archive and adding them to reprocessing log.\n')
                sys.stdout.flush()
                qalas_in_qc_but_not_archive.append(temp_session)
                continue


            sys.stdout.write('   QALAS files found within: {}\n'.format(qalas_folders))
            sys.stdout.flush()

            sys.stdout.write('   Supplemental info associated with archive: {}\n'.format(supplemental_infos))
            sys.stdout.flush()
            
            #Run SyMRI to create relaxometry maps and convert to BIDS
            output_info = convert_single_tar(qalas_folders, supplemental_infos, best_qalas_info,
                                            session_base_dir,
                                            args.symri_container_path,
                                            symri_layout,
                                            symri_global,
                                            dcm2bids_config)
            
            if output_info['num_niftis_generated'] < 3:
                no_niftis_to_upload.append(temp_session)
                sys.stdout.write('   Participant didnt have any data to upload... saving info for debugging.')
                sys.stdout.flush()
                continue
            
            sys.stdout.write('   Info associated with unpacking: {}\n'.format(output_info))
            sys.stdout.flush()
            
            #Send the results to s3
            status = push_to_s3(output_info['bids_path'], output_info['subject_label'], bucket_name = args.custom_loris_bucket_name,
                                    prefix = os.path.join('derivatives', 'ses-' + output_info['session_label'], 'symri'),
                                    different_config_path=args.loris_config_path)    

            if status == False:
                raise ValueError('Error: Pushing data to S3 was unsuccessful.')
            
            sys.stdout.write('   Data successfully uploaded to S3: {}\n'.format(output_info))
            sys.stdout.flush()
            
            #After BIDS conversion + uploading to S3, update the
            #tracking dictionary with data for this subject.
            tracking_log_data[temp_session] = best_qalas_info

            #If this session was in the reprocess list, remove it
            if temp_session in reprocess_attempted:
                reprocess_attempted.remove(temp_session)
                if temp_session in new_archive_available: new_archive_available.remove(temp_session)
                if temp_session in qalas_in_qc_but_not_archive: qalas_in_qc_but_not_archive.remove(temp_session)
                if temp_session in no_niftis_to_upload: no_niftis_to_upload.remove(temp_session)
                if temp_session in missing_archive: missing_archive.remove(temp_session)
        
        else:
            
            #Add the jsons for this session to the json log so they
            #can be ignored in further processing.
            if qc_or_other_missing == False:
                best_qalas_info = {}
                best_qalas_info['jsons_tested'] = jsons_dict[temp_session]
                tracking_log_data[temp_session] = best_qalas_info
                sys.stdout.write('   No signs of QALAS scan that should be used for processing. Skipping subject archive and adding them to log.\n')
                sys.stdout.flush()
            else:
                sys.stdout.write('   Processing of archive has been skipped due to missing QC information. Subject info will not be added to log.\n')
                sys.stdout.flush()

        #Clean up the session working directory
        if args.keep_work_dirs == False:
            os.chdir(args.base_directory_for_proc)
            shutil.rmtree(session_base_dir)
            sys.stdout.write('   Removing working directory at: {}\n\n\n'.format(session_base_dir))
            sys.stdout.flush()
            
    #Save the tracking log
    date_time = datetime.now().strftime("date%Y_%m_%d_time%H_%M_%S")
    new_path = os.path.join(logs_dir, 'tracking_log_{}.json'.format(date_time))
    with open(new_path, 'w') as f:
        f.write(json.dumps(tracking_log_data, indent = 5))
    sys.stdout.write('Tracking log has been updated at: {}\n'.format(new_path))
    sys.stdout.flush()


    #Save the reprocess log
    reproc_log_data['to_reprocess'] = [i for i in to_reprocess if i not in reprocess_attempted]
    new_archive_available = list(set(new_archive_available))
    reproc_log_data['new_archive_available'] = new_archive_available
    reproc_log_data['qalas_in_qc_but_not_archive'] = list(set(qalas_in_qc_but_not_archive))
    reproc_log_data['no_niftis_to_upload'] = list(set(no_niftis_to_upload))
    reproc_log_data['missing_archive'] = list(set(missing_archive))

    with open(reproc_log_path, 'w') as f:
        f.write(json.dumps(reproc_log_data, indent = 5))
    sys.stdout.write('Reprocessing log has been updated at: {}\n'.format(reproc_log_path))
    sys.stdout.flush()


    sys.stdout.write('\n\nCurrent batch of processing has been completed. Exiting processing now.\n')
    sys.stdout.write('Reprocessing was attempted for: {}'.format(reprocess_attempted))
    sys.stdout.flush()

    return


if __name__ == "__main__":
    main()


