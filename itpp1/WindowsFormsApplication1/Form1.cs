using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace WindowsFormsApplication1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            var start = DateTime.Now;
            var p = YourWrapperClass.search("C:\\Users\\yue\\Documents\\Visual Studio 2015\\Projects\\itpp1\\Win32\\Debug\\f1815.3_s19.2_bw20_0.08s_hackrf-1.bin");
            string res = Marshal.PtrToStringAnsi(p);
            var end = DateTime.Now;
            Console.WriteLine("Result from C++ method: " + res);
            Console.WriteLine("time: " + (end - start).TotalSeconds);
            Console.WriteLine("Result from C++ method: " + res);

            // YourWrapperClass.search();
            // Console.WriteLine("Result from C++ method: ");
        }
    }
}
