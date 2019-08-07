create_category_database = __import__('1_create_category_database')
edit_raw_XML = __import__('2_edit_raw_XML')
parse_to_dlookup = __import__('3_parse_to_dlookup')
parse_to_qlookup = __import__('4_parse_to_qlookup')
parse_to_tlookup = __import__('5_parse_to_tlookup')
parse_to_jlookup = __import__('6_parse_to_jlookup')
create_tutorial_coll = __import__('7_create_tutorial_coll')
import helpers as hp

def main():
    root = hp.get_path()
    hp.sign_off()
    hp.space()

    ## This script will only work with python 3. ##

    ## Stage 0
    ## There is a prior stage here where make_object_list.py is called by Xcode.
    ## This produces the header file which Build_Max_Docs.cpp uses to know about FrameLib objects and types.
    
    ## Creates a category database in .json format.
    ## The JSON file is used by edit_raw_XML.py to assign object categories to the xml files.
    print('1. Building Category Database')
    create_category_database.main(root)
    hp.hyp()

    ## The purpose of this script is to set the categories for the Raw_XML files. 
    ## C++ doesnt know about the categories at XML creation and its easier to iterate file structures in python.
    ## Edited XML files are copied from /tmp/ to the refpages directory
    print('2. Editing XML Files')
    edit_raw_XML.main(root)
    hp.hyp()

    ## This script creates a dictionary used to display specific object info in the extras Max Patch.
    ## Similar to the qlookup, but is specifically used to display the digest with mouse hovering
    print('3. Building dlookup')
    parse_to_dlookup.main(root)
    hp.hyp()

    ## This script creates a dictionary that contains specific object information.
    ## This provides the dynamic hover behaviour
    print('4. Building qlookup')
    parse_to_qlookup.main(root)
    hp.hyp()

    ## Creates a dictionary used to display names and descriptions of tutorials in the extras Max Patch.
    ## The tutorials are categorised by difficulty. {Beginner, Intermediate, Advanced}
    print('5. Building tlookup')
    parse_to_tlookup.main(root)
    hp.hyp()

    ## Creates a dict containing information about object parameters. This is used by the help file template.
    print('6. Building jlookup')
    parse_to_jlookup.main(root)
    hp.hyp()

    ## Creates a coll containing the file names of the tutorials. Makes it a bit easier to load them.
    print('7. Building tutorial name coll')
    create_tutorial_coll.main(root)
    hp.hyp()
    print(' ')
    print("Completed all python scripts!!!")


if __name__ == '__main__':
    main()

