"""
---
    Copyright (c) 2018 Baskar Ganapathysubramanian, Balaji Sesha Sarath Pokuri
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
---
"""

## --- end license text --- ##
import tarfile
from distutils.version import LooseVersion
import re
import tempfile
import shlex
import os
from functools import partial
from typing import List, Union
from stat import S_ISDIR

import paramiko
import socket

# use easygui if it's available, so we don't have passwords echoed when using PyCharm
try:
    from easygui import passwordbox
    GET_PRIVATE_INPUT = passwordbox
except ImportError:
    from getpass import getpass
    GET_PRIVATE_INPUT = getpass

import sys

# use dill for call_on_remote, which is optional
try:
    from dill import dill
except ImportError:
    dill = None


class Host:
    def __init__(self, username, hostname, port=22):
        self.username = username
        self.hostname = hostname
        self.port = port
        self.agent = paramiko.Agent()

    def get_password(self):
        """

        :return: the password for the host
        """
        return GET_PRIVATE_INPUT('Password: ')

    def get_keys(self):
        """

        :return: a list of public/private keys to try authenticating with
        """
        return self.agent.get_keys()

    def get_interactive(self, title: str, instructions: str, prompts: List[str]) -> List[str]:
        """
        Handles the ssh 'interactive' authentication mode (user answers a series of prompts).

        :param title: title of the window
        :param instructions: instructions, to be shown before any prompts
        :param prompts: the list of prompts
        :return: the list of responses (in the same order as the prompts)
        """
        if instructions:
            print(instructions)

        results = []
        for prompt in prompts:
            results.append(GET_PRIVATE_INPUT(prompt[0]))
            # if prompt[1]:
            #     results.append(GET_PRIVATE_INPUT(prompt[0]))
            # else:
            #     results.append(input(prompt[0]))

        return results


class Connection:
    def __init__(self):
        self.host = None  # type: Host
        self._transport = None  # type: paramiko.Transport
        self._sftp = None  # type: paramiko.SFTPClient
        self._remote_python = None  # type: str

    def connect(self, host: Host):
        self.host = host
        self._transport = self._open_connection()
        self._sftp = self._transport.open_sftp_client()

    def sftp(self) -> paramiko.SFTPClient:
        return self._sftp

    def put_file(self, local_path: str, remote_path: str) -> None:
        """
        Uploads a local file to the remote host, same as paramiko.SFTPClient.put
        """
        self._sftp.put(local_path, remote_path)

    def mkdirs(self, remote_dir) -> None:
        """
        Creates remote_dir recursively as a directory on the remote, creating any necessaries parent directories
        along the way. Similar to `mkdir -p`, except it doesn't error if the directory already exists.
        """
        sftp = self.sftp()
        dirs_ = []
        dir_ = remote_dir
        while len(dir_) > 1:
            dirs_.append(dir_)
            dir_, _ = os.path.split(dir_)

        if len(dir_) == 1 and not dir_.startswith("/"):
            dirs_.append(dir_)  # For a remote path like y/x.txt

        while len(dirs_):
            dir_ = dirs_.pop()
            try:
                sftp.stat(dir_)
            except IOError:
                sftp.mkdir(dir_)

    def put_dir(self, local_path: str, remote_path: str) -> None:
        """
        Compresses local_path into a .tar.gz archive, uploads it to the remote, extracts it into remote_path,
        and finally deletes the temporary tar archive. Assumes the remote has the 'tar' utility available.
        """

        # first, create the folder on the remote
        self.mkdirs(remote_path)

        with tempfile.TemporaryFile() as f:
            # write tar file
            with tarfile.open(fileobj=f, mode='w:gz') as tarf:
                for root, dirs, files in os.walk(local_path):
                    for file in files:
                        p = os.path.join(root, file)
                        rel_p = os.path.relpath(p, local_path)
                        tarf.add(p, arcname=rel_p)

            # transfer to remote
            remote_archive_path = remote_path.rstrip('/') + '_put.tar.gz'
            f.seek(0)  # move read cursor to start of file
            self._sftp.putfo(f, remote_archive_path)

        # unzip on remote and remote the zip file
        cmd = 'tar xf {} --directory {} && rm {}'.format(
            shlex.quote(remote_archive_path), shlex.quote(remote_path), shlex.quote(remote_archive_path))
        self.exec_command(cmd)

    # from www.stackoverflow.com/questions/24427283/getting-a-files-from-remote-path-to-local-dir-using-sftp-in-python
    def get_dir(self, remote_dir: str, local_dir: str) -> None:
        """
        Download directory from remote directory to local directory
        """
        dir_items = self.sftp().listdir_attr(remote_dir)
        for item in dir_items:
            remote_path = os.path.join(remote_dir, item.filename)
            local_path = os.path.join(local_dir, item.filename)
            if S_ISDIR(item.st_mode):
                self.get_dir(remote_path, local_path)
            else:
                self.sftp().get(remote_path, local_path)

    def get_file(self, remote_path: str, local_path: str) -> None:
        """
        Downloads a file from the remote, same as paramiko.SFTPClient.get
        """
        self._sftp.get(remote_path, local_path)

    def exec_command(self, cmd: str, cwd=None, check_exitcode=True, encoding='utf-8') -> (str, str, int):
        """
        Executes cmd in a new shell session on the remote host.

        :param cmd: command to execute
        :param cwd: directory to execute the command in - performed by prepending 'cd [cwd] && ' to cmd
        :param check_exitcode: if true, instead of returning the exit code of cmd as part of the return tuple, verify
                that the return code is zero. If it is not, an exception is raised with the contents of stderr.
        :param encoding: encoding to decode stdout/stderr with. Defaults to utf-8.
        :return: if check_exitcode is True, (stdout: str, stderr: str). If it is False, (stdout, stderr, rc: int).
                stdout and stderr are decoded according to encoding.
        """
        tp = self._transport
        assert tp is not None

        if cwd is not None:
            cmd = 'cd ' + shlex.quote(cwd) + ' && ' + cmd

        channel = tp.open_session()  # type: paramiko.Channel
        stdoutf = channel.makefile('r')
        stderrf = channel.makefile_stderr('r')
        channel.exec_command(cmd)

        stdout = stdoutf.read().decode(encoding)
        stderr = stderrf.read().decode(encoding)
        exitcode = channel.recv_exit_status()
        stdoutf.close()
        stderrf.close()
        channel.close()

        if check_exitcode and exitcode != 0:
            raise Exception("Bad exit status code ({}) - {}".format(exitcode, stderr))

        if not check_exitcode:
            return stdout, stderr, exitcode
        else:
            return stdout, stderr

    def remote_python(self) -> str:
        """
        Returns a string that, when invoked as a command on the remote, will execute a Python that:

        * Matches the version that this script was invoked with (i.e. matching sys.version_info)
        * Has the 'dill' module installed

        The remote Python is discovered by trial and error using common Python names.
        The search is performed once and then cached.
        If no such Python is available, this will return None.
        """
        if not self._remote_python:
            self._remote_python = self._detect_remote_python()
        return self._remote_python

    def call_on_remote(self, remote_func, *args, remote_cwd: Union[str, None]=None):
        """
        Call a function created on this system on a remote system with args.
        This is done by pickling it with dill, SFTPing it to a file on the remote, executing a Python script on the
        remote that un-dills the file, calls the function, dills the result and prints it to stdout.
        Finally, stdout is un-dilled on the local machine to give the return value.
        This requires the remote to have a matching Python version.

        :param remote_func: function to call
        :param args: any arguments to call the function with
        :param remote_cwd: directory on the remote to call the script from (must have write access to this directory)
        :return: value returned by f
        """
        if not dill:
            raise Exception('dill not installed - cannot use call_on_remote.')

        if self.remote_python() is None:
            raise Exception('Remote Python 3 installation not found - cannot use call_on_remote.')

        # serialize remote_func using dill and save it to a file on the remote
        remote_func = partial(remote_func, *args)
        if remote_cwd is not None:
            remote_path = os.path.join(remote_cwd, 'exec.dill')
        else:
            remote_path = 'exec.dill'

        with tempfile.TemporaryFile() as fl:
            dill.dump(remote_func, fl)
            fl.seek(0)  # move read cursor to start of file
            self._sftp.putfo(fl, remote_path)

        # load remote_func using dill on the remote and call it
        cmd = self.remote_python() \
              + " -c 'from dill import dill; import sys; sys.stdout.buffer.write(dill.dumps(dill.load(open(\"exec.dill\", \"rb\"))()))'"
        if remote_cwd is not None:
            cmd = 'cd ' + shlex.quote(remote_cwd) + ' && ' + cmd

        tp = self._transport
        channel = tp.open_session()  # type: paramiko.Channel
        stdoutf = channel.makefile('rb')
        stderrf = channel.makefile_stderr('r')
        channel.exec_command(cmd)
        exitcode = channel.recv_exit_status()

        if exitcode == 0:
            ret = dill.load(stdoutf)
            stdoutf.close()
            stderrf.close()
            channel.close()
            return ret
        else:
            err = stderrf.readlines()
            stdoutf.close()
            stderrf.close()
            channel.close()
            raise Exception("Remote execution error:\n\t" + "\t".join(err))

    def _open_connection(self) -> paramiko.Transport:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((self.host.hostname, self.host.port))
        tp = paramiko.Transport(sock)
        tp.set_keepalive(120)  # keepalive every 2m
        tp.auth_timeout = 150  # time until authentication will wait before timeout -- default : 30 sec
        tp.default_max_packet_size = 100000000
        tp.default_window_size = 10000000
        tp.start_client()

        # Authenticate
        # first try auth_none so we have a list of auths to check
        auth_types = []
        try:
            tp.auth_none(self.host.username)
            return tp  # that somehow worked
        except paramiko.BadAuthenticationType as e:
            auth_types = e.allowed_types

        # try the auth types that the server said are available
        for auth_type in auth_types:
            try:
                if auth_type == 'publickey':
                    keys = self.host.get_keys()
                    if len(keys) == 0:
                        continue

                    accepted = False
                    for key in keys:
                        try:
                            tp.auth_publickey(self.host.username, key)
                            accepted = True
                            break
                        except paramiko.AuthenticationException:
                            pass
                    if not accepted:
                        raise paramiko.AuthenticationException()

                elif auth_type == 'password':
                    tp.auth_password(self.host.username, self.host.get_password())
                elif auth_type == 'interactive':
                    tp.auth_interactive(self.host.username, self.host.get_interactive)
                elif auth_type == 'keyboard-interactive':
                        tp.auth_interactive(self.host.username, self.host.get_interactive)
                else:
                    print("Skipping authentication type '" + auth_type + "'")

            except paramiko.AuthenticationException as e:
                print("Authentication type '" + auth_type + "' failed")
                print("Exception: --" + e.__str__())

            # authenticated successfully
            if tp.is_authenticated():
                return tp

        raise paramiko.AuthenticationException()

    def _detect_remote_python(self,
                              req_ver: LooseVersion=None,
                              required_modules: List[str] = list(['dill'])) -> Union[str, None]:
        if req_ver is None:
            vs = [sys.version_info.major, sys.version_info.minor, sys.version_info.micro]
            vs = [str(s) for s in vs]
            req_ver = LooseVersion('.'.join(vs))

        guesses = ['python', 'python3', 'python3.5', 'python3.6',
                   'module load python && python3.5', 'module load python && python']
        for guess in guesses:
            print("Testing for remote Python '" + guess + "'")
            stdout, stderr, rc = self.exec_command(guess + ' --version', check_exitcode=False)
            if rc != 0:
                print("  Not found")
                continue

            output = stdout if len(stdout) > 0 else stderr
            match = re.match(r"Python (\d+\.\d+\.\d+)", output)
            if not match:
                print("  Invalid version format (stdout: " + stdout + ', stderr: ' + stderr + ")")
                continue

            ver = LooseVersion(match.group(1))
            if ver != req_ver:
                print('  Version ' + str(ver) + ' does not match required (' + str(req_ver) + ')')
                continue

            ok = True
            for modname in required_modules:
                if not self._check_remote_python_module(guess, modname):
                    print("  Missing required module '" + modname + "'")
                    ok = False
                    break

            if not ok:
                continue

            print('  OK')
            return guess

        return None

    def _check_remote_python_module(self, python: str, modname: str):
        """
        check if module 'modname' is installed on the remote using the given python command
        """
        _, _, rc = self.exec_command(python + ' -c "import ' + modname + '"', check_exitcode=False)
        return rc == 0
