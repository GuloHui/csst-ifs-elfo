import setuptools

# 读取README.md作为长描述
with open("README.md", "r") as f:
    long_description = f.read()

# 读取依赖列表requirements.txt
# 忽略#开头或者版本号不明确指定的条目
with open("requirements.txt", "r") as f:
    requirements = [
        req.strip()
        for req in f.readlines()
        if not req.startswith("#") and req.__contains__("==")
    ]

# 配置、安装
setuptools.setup(
    name="elfo",  # 包名
    version="0.1.0",  # 版本号
    author="Hui Guo",  # 作者
    author_email="guohui@mail.ustc.edu.cn",  # 邮箱
    description="The CSST pipeline - IFS elfo",  # 短描述
    long_description=long_description,  # 长描述
    long_description_content_type="text/markdown",  # 长描述类型
    url="https://github.com/GuloHui/csst-ifs-elfo",  # 主页
    packages=["src/csst_ifs_elfo"],  # 用setuptools工具自动发现带有__init__.py的包
    license="GNU General Public License v3.0",  # 证书类型
    classifiers=[  # 程序分类, 参考 https://pypi.org/classifiers/
        # How mature is this project?
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    include_package_data=True,  # 设置包含随包数据
    package_data={  # 具体随包数据路径
        },
    # 请注意检查，防止临时文件或其他不必要的文件被提交到仓库，否则会一同安装
    python_requires=">=3.11",  # Python版本要求
    install_requires=requirements,
)