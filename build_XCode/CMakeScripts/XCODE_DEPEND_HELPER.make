# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.shellbenchmark_bin.Debug:
PostBuild.igl.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.igl_opengl_glfw.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.igl_opengl_glfw_imgui.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.alglib.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.igl_opengl_glfw.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.igl_opengl.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.igl.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.igl_common.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.imgui.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.glfw.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
PostBuild.glad.Debug: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin:\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Debug/libalglib.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Debug/libimgui.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Debug/libglfw3.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Debug/libglad.a
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Debug/shellbenchmark_bin


PostBuild.glad.Debug:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Debug/libglad.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Debug/libglad.a


PostBuild.glfw.Debug:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Debug/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Debug/libglfw3.a


PostBuild.imgui.Debug:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Debug/libimgui.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Debug/libimgui.a


PostBuild.alglib.Debug:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Debug/libalglib.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Debug/libalglib.a


PostBuild.shellbenchmark_bin.Release:
PostBuild.igl.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.igl_opengl_glfw.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.igl_opengl_glfw_imgui.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.alglib.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.igl_opengl_glfw.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.igl_opengl.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.igl.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.igl_common.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.imgui.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.glfw.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
PostBuild.glad.Release: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin:\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Release/libalglib.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Release/libimgui.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Release/libglfw3.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Release/libglad.a
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/Release/shellbenchmark_bin


PostBuild.glad.Release:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Release/libglad.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Release/libglad.a


PostBuild.glfw.Release:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Release/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Release/libglfw3.a


PostBuild.imgui.Release:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Release/libimgui.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Release/libimgui.a


PostBuild.alglib.Release:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Release/libalglib.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Release/libalglib.a


PostBuild.shellbenchmark_bin.MinSizeRel:
PostBuild.igl.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.igl_opengl_glfw.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.igl_opengl_glfw_imgui.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.alglib.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.igl_opengl_glfw.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.igl_opengl.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.igl.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.igl_common.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.imgui.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.glfw.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
PostBuild.glad.MinSizeRel: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin:\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/MinSizeRel/libalglib.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/MinSizeRel/libimgui.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/MinSizeRel/libglfw3.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/MinSizeRel/libglad.a
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/MinSizeRel/shellbenchmark_bin


PostBuild.glad.MinSizeRel:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/MinSizeRel/libglad.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/MinSizeRel/libglad.a


PostBuild.glfw.MinSizeRel:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/MinSizeRel/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/MinSizeRel/libglfw3.a


PostBuild.imgui.MinSizeRel:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/MinSizeRel/libimgui.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/MinSizeRel/libimgui.a


PostBuild.alglib.MinSizeRel:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/MinSizeRel/libalglib.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/MinSizeRel/libalglib.a


PostBuild.shellbenchmark_bin.RelWithDebInfo:
PostBuild.igl.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.igl_opengl_glfw.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.igl_opengl_glfw_imgui.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.alglib.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.igl_opengl_glfw.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.igl_opengl.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.igl.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.igl_common.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.imgui.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.glfw.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
PostBuild.glad.RelWithDebInfo: /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin:\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/RelWithDebInfo/libalglib.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/RelWithDebInfo/libimgui.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/RelWithDebInfo/libglfw3.a\
	/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/RelWithDebInfo/libglad.a
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/RelWithDebInfo/shellbenchmark_bin


PostBuild.glad.RelWithDebInfo:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/RelWithDebInfo/libglad.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/RelWithDebInfo/libglad.a


PostBuild.glfw.RelWithDebInfo:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/RelWithDebInfo/libglfw3.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/RelWithDebInfo/libglfw3.a


PostBuild.imgui.RelWithDebInfo:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/RelWithDebInfo/libimgui.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/RelWithDebInfo/libimgui.a


PostBuild.alglib.RelWithDebInfo:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/RelWithDebInfo/libalglib.a:
	/bin/rm -f /Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/RelWithDebInfo/libalglib.a




# For each target create a dummy ruleso the target does not have to exist
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Debug/libalglib.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/MinSizeRel/libalglib.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/RelWithDebInfo/libalglib.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/alglib/Release/libalglib.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Debug/libglad.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/MinSizeRel/libglad.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/RelWithDebInfo/libglad.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glad/Release/libglad.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Debug/libglfw3.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/MinSizeRel/libglfw3.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/RelWithDebInfo/libglfw3.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/glfw/src/Release/libglfw3.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Debug/libimgui.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/MinSizeRel/libimgui.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/RelWithDebInfo/libimgui.a:
/Users/chenzhen/UT/Research/Projects/shellbenchmarks/build_XCode/imgui/Release/libimgui.a:
