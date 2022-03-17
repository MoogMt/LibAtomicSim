// Guard
#ifndef generalio_h
#define generalio_h

using namespace std;

bool fileExistsIO( ifstream &file_handler );
bool fileExists( string file_path );

int getNumberLineIO( ifstream &file_handler );
int getNumberLine( string file_path );

#endif