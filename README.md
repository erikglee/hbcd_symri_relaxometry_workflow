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
        consider_reprocessing portion of the json, any new archives for the
        session will be ignored. If the user has added a subject/session
        to the to_reprocess portion of the json, then the code will attempt
        to overwrite any existing processing for that subject/session. This
        may result in either the same or different processing results depending
        on the new data that is available for the subject. If the reproc_log.json
        file does not exist, then the code will create it.
    3. Iterate through all subjects/sessions that need to be processed. For subjects/
        sessions that have a QALAS scan and are ready to be processed, download the
        archive containing the best QALAS scan, unpack the archive, and identify
        the dicom folder that contains the QALAS scan.
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