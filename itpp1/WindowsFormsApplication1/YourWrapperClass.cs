using System;
using System.Runtime.InteropServices;
using System.Text;

namespace WindowsFormsApplication1
{
    public class YourWrapperClass
    {
        [DllImport("C:\\Users\\yue\\Documents\\Visual Studio 2015\\Projects\\itpp1\\Win32\\Debug\\itpp1.dll", EntryPoint = "search", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr search(string a);
        // public static extern void search();
    }
}
