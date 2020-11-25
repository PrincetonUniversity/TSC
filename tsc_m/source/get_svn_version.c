const char* get_svn_version_(void) 
{ const char* SVN_Version = "245M";
  int my_svn_number = atoi("245M");
  static int count = 0; if(!count++) printf(" TSC Version : %s\n", SVN_Version);
  return SVN_Version; 
}
