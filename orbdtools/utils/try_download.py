import wget

def wget_download(url,dir_file,desc=None):
    """
    Download files by wget command

    Inputs:
        url -> [str] URL of the file to be downloaded
        dir_file -> [str] Path of the file to store
        desc -> [str,optional,default=None] Description of the downloading   
    Outputs:
        wget_out -> [str] Path of the file downloaded

    """
    if desc: print(desc)
    wget_out = wget.download(url,dir_file)
    print()

    return wget_out