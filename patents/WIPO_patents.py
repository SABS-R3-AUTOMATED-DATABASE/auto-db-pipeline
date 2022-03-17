import os
import ftplib
import zipfile


def get_seq_listing(URL: str):
    """
    Access the ftp server of WIPO patent to download Sequence listing files uploaded with the patent
    Download and read the txt files and add the sequences into Content column
    """
    output = ""
    year = URL[36:40]
    wo_folder = URL[34:36] + URL[38:40] + "_" + URL[40:46]
    ftp = ftplib.FTP("ftp.wipo.int")
    ftp.login()
    path = "/pub/published_pct_sequences/publication/"
    ftp.cwd(path + year + "/")
    filelist = [item for item in ftp.nlst() if "." not in item]
    for file in filelist:
        ftp.cwd(path + year + "/" + file + "/")
        filelist2 = [item for item in ftp.nlst() if "." not in item]
        if wo_folder in filelist2:
            ftp.cwd(wo_folder)
            ziplist = [file for file in ftp.nlst() if file != "applicant.txt"]
            if not os.path.exists("/data/temp"):
                os.makedirs("/data/temp")
            else:
                for item in os.listdir("/data/temp"):
                    os.remove("/data/temp/" + item)
            for zip in ziplist:
                ftp.retrbinary("RETR " + zip, open("/data/temp/" + zip, "wb").write)
                with zipfile.ZipFile("/data/temp/" + zip, "r") as zip_ref:
                    zip_ref.extractall("/data/temp/")
            ftp.close()
            txtlist = [file for file in os.listdir("/data/temp") if ".txt" in file]
            for txt in txtlist:
                with open("/data/temp/" + txt, "r", errors="ignore") as f:
                    sl = [line.rstrip("\n") for line in f]
                    output = output + " ".join(sl)
            for item in os.listdir("/data/temp"):
                os.remove("/data/temp/" + item)
            os.rmdir("/data/temp")
            break
    return output


def get_wipo_sequences(df):
    """
    Check if there is potential sequence in WIPO patents, if so, download possible sequence listing file and add to Content column
    """
    for i in range(df.shape[0]):
        if "WO" in df.loc[i, "URL"]:
            dummy = False
            for _ in df.loc[i, "Claim"]:
                if "seq id no" in _.lower():
                    dummy = True
                    break
            if dummy:
                seq_list = get_seq_listing(df.loc[i, "URL"])
                df.loc[i, "Content"] = df.loc[i, "Content"] + [seq_list]
    return df
