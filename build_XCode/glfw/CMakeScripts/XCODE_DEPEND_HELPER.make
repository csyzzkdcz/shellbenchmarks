# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.glfw.Debug:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/Debug/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/Debug/libglfw3.a


PostBuild.glfw.Release:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/Release/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/Release/libglfw3.a


PostBuild.glfw.MinSizeRel:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/MinSizeRel/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/MinSizeRel/libglfw3.a


PostBuild.glfw.RelWithDebInfo:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/RelWithDebInfo/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_Xcode/glfw/src/RelWithDebInfo/libglfw3.a




# For each target create a dummy ruleso the target does not have to exist
