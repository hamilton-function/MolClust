from PIL import Image, ImageDraw, ImageFont
import csv
import pubchempy as pcp
import os


'''
reads in a csv of structure:
col0: cid 
col1: compound name

returns: dictionary with cid as key and compound name as value
'''
def get_cid_compound_dict(csvfile):
    ciddict = {}
    with open(csvfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            ciddict[row[0]] = row[1]
    return ciddict



'''
routine which takes a png file and applies text (compound name) on it and saves it as a new file
'''
def compname2image(compname,pic):
    # create Image object with the input image
    image = Image.open(pic)

    # initialise the drawing context with the image object as background
    draw = ImageDraw.Draw(image)

    # create font object with the font file and specify desired size
    font = ImageFont.truetype("C:\\WINDOWS\\FONTS\\MICROSS.TTF", 16)

    # draw the message on the background
    draw.text((0, 0), compname, fill=150, font=font)

    #imagename
    path, file = os.path.split(pic)

    # save the edited image
    image.save(path+'\\'+file)




'''
now do the naming for all images
'''
def compoundsImagesWithNames(ciddict, targetdir):
    for cid in ciddict:
        imgpath = targetdir + str(cid) + '.png'
        pcp.download('PNG', targetdir + str(cid) + '.png', str(cid), 'cid', overwrite=True)
        compname = ciddict[cid]
        compname2image(compname, imgpath)


cidict = get_cid_compound_dict('E:\\MolClust\\covidmancandidates.csv')
compoundsImagesWithNames(cidict, "E:\\MolClust\\namedcompoundimages\\")