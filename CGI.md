
#STRINGMOL WEB API

*From Adam Nellis:*

We need three entry points into the Stringmol code:

1. bind: A function that determines if/where two molecules bind.
- step: A function that runs one step of a reaction.
- bucket: A function that starts a bucket running.

I had these implemented as three separate CGI scripts that had their inputs encoded as GET parameters in the url (like this: /path/script_name.cgi?param1=value1&param2=value2). 
Their output was a single string, encoded in the same way (output1=value1&output2=value2). This wasn't ideal, and should be improved. For example, could return the output as a json string.

Below, I've written down the inputs and outputs of these three functions. 
This is how I implemented them, but there is some redundancy in there that we could get rid of. Feel free to improve them if you want.

##bind
###Performs:
Work out if and where two molecules bind.

###Inputs:
    String - molecule A
    String - molecule B
 
###Outputs:
    Boolean - Whether the bind is possible
    Int - Length of the bind
    Int - Probability of the bind
    Int - Position of the start of the bind on molecule A
    Int - Position of the end of the bind on molecule A
    Int - Position of the start of the bind on molecule B
    Int - Position of the end of the bind on molecule B
    String - molecule A
    String - molecule B
    Int - Position of molecule A's instruction pointer
    Int - Position of molecule B's instruction pointer
    Int - Position of molecule A's flow pointer
    Int - Position of molecule B's flow pointer
    Int - Position of molecule A's read pointer
    Int - Position of molecule B's read pointer
    Int - Position of molecule A's write pointer
    Int - Position of molecule B's write pointer
    Boolean - if the instruction pointer is active on molecule A
    Boolean - if the flow pointer is active on molecule A
    Boolean - if the read pointer is active on molecule A
    Boolean - if the write pointer is active on molecule A


##Step

###Performs:
    Run one step of an existing bind.

###Inputs: 
    String - molecule A
    String - molecule B
    Int - position of molecule A's instruction pointer
    Int - position of molecule B's instruction pointer
    Int - position of molecule A's flow pointer
    Int - position of molecule B's flow pointer
    Int - position of molecule A's read pointer
    Int - position of molecule B's read pointer
    Int - position of molecule A's write pointer
    Int - position of molecule B's write pointer
    Boolean - if the instruction pointer is active on molecule A
    Boolean - if the flow pointer is active on molecule A
    Boolean - if the read pointer is active on molecule A
    Boolean - if the write pointer is active on molecule A
 
###Outputs:
    String - molecule A
    String - molecule B
    Int - position of molecule A's instruction pointer
    Int - position of molecule B's instruction pointer
    Int - position of molecule A's flow pointer
    Int - position of molecule B's flow pointer
    Int - position of molecule A's read pointer
    Int - position of molecule B's read pointer
    Int - position of molecule A's write pointer
    Int - position of molecule B's write pointer
    Boolean - if the instruction pointer is active on molecule A
    Boolean - if the flow pointer is active on molecule A
    Boolean - if the read pointer is active on molecule A
    Boolean - if the write pointer is active on molecule A
    String - molecule that has been cleaved off at this step
    Boolean - if the reaction has ended


##Bucket

###Performs:
    
Set a bucket running.

###Inputs:
    Int - run number (unique ID)
    Bag - molecules in the bucket
    Float - cell radius
    Float - agent radius
    Int - initial energy
    Boolean - mutation on or off
    
The bag of molecules is represented as a string formatted like:
    "(%d,%s)(%d,%s)(%d,%s)..."
    so each molecule species is a pair: (number of copies, molecule string)

###Outputs:
    
	None.
Stringmol is started asynchronously, and this function returns straight away.
    Stringmol writes its results to a directory named by the run number.
    A different script (in PHP) polls the results directory for this run.
    The stringmol run can be stopped early, by writing a file called "_STOP" in the results directory (Stringmol checks for this file existing every few iterations).



##Notes

The way Adam did it was to have the CGI script write to stdout. When JavaScript calls it, the response it gets back is whatever the script wrote to stdout. I didn't encode it as HTML. I encoded it in the same format as the input string: param1=value1&param2=value2...
But I didn't like that, so we could improve it. I would do the output as json, because JavaScript can parse that easily, but I'm happy with whatever you think best.

Also, I don't know how much development you want to do on this. The simplest way would be to create three CGI scripts like I did before. But if we want to be more ambitious, then we could make a Stringmol library that exposes these three functions. Then we could have PHP scripts that use the library. This would make it easier to run Stringmol from other languages, and easier for other people to use Stringmol in their own applications.

Regards,
Adam
